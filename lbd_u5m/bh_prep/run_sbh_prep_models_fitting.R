## ##############################################################################
## ##############################################################################
## PURPOSE: Run SBH Model for selected region
##           Launched from a qsub with region and run_date as args
##           Model fit objects get used later for prediction of SBH-only data
## INPUTS:  region argument, ab_times, CBH/SBH dataset
##
## ##############################################################################
## ##############################################################################



## SGE LOG HEADER
message(' ... ... ... ... ... ... ... ... ... ... ')
message('Run notes:')
print(R.Version())
print(Sys.info())
message(' ... ... ... ... ... ... ... ... ... ... ')


 
##############################################################################################
message('\nSETUP ... ')

# Load in args from qsub call
run_date   <- as.character(commandArgs()[4])
reg        <- as.character(commandArgs()[5])

# print out message describing arg inputs for this model run
message(sprintf('\nRUN_DATE: %s',run_date))
message(sprintf('\nREGION: %s',reg))

# Get user, assume code is in users directory
user      <- Sys.info()['user']
message(sprintf('\nUSER: %s',user))

# Get output directory
root   <- '<<<< FILEPATH REDACTED >>>>'

core_repo <- sprintf("<<<< FILEPATH REDACTED >>>>", user)
source(paste0(core_repo, '<<<< FILEPATH REDACTED >>>>'))

# Load packages
packlist <- c('data.table','parallel','tictoc','seegMBG','survival','magrittr','scales','lme4','Biograph')
load_R_packages(packlist)

# load custom functions
setwd(sprintf('<<<< FILEPATH REDACTED >>>>', user))
source('sbh_utils.R')
message(' ... Custom functions sucessfully loaded')

# make directories
make_output_dirs(root = root, run_date = run_date, country = reg)
rd_dir   <- sprintf('%s/%s', root, run_date)
rdc_dir  <- sprintf('%s/%s', rd_dir, reg)
dat_dir  <- '<<<< FILEPATH REDACTED >>>>'

## load configs, which were saved in launch script
# ab times
ab_times <- fread(sprintf('<<<< FILEPATH REDACTED >>>>',rd_dir))

# config values into environment
load_config(rd_dir)

##############################################################################################


##############################################################################################
message('\nLOAD IN CBH-SBH COMBINED DATASET ... ')
tic.clearlog()
tic('data_loading')

# load in the massive DHS dataset
data <- setDT(readRDS("<<<< FILEPATH REDACTED >>>>"))
info <- fread("<<<< FILEPATH REDACTED >>>>")
data <- merge(data, info, by = 'country', all.x = TRUE)

# map between region names
regmap <- data.table(MAP_gbdregion = c("North Africa / Middle East", "Sub-Saharan Africa, West/Central", "Asia", 
                                       "Sub-Saharan Africa, South/East", "Latin America and the Caribbean"),
                     region        = c('NAME','SSAWC','ASIA','SSASE','LAC'))
data <- merge(data,regmap,by='MAP_gbdregion',all.x=TRUE)

# subset out the DHS/MIS data only
data <- subset(data, survey %in% c('MACRO_DHS','MACRO_MIS','NIC/DHS_ENDESA'))
data <- subset(data,! nid %in% 20786)
data$nid <- as.factor(data$nid)

message(sprintf(' ... Full database with %i rows sucessfully loaded ',nrow(data)))
toc(log = TRUE) # end data_loading tic

##############################################################################################


##############################################################################################
message('\nLOAD SDI INFORMATION ... ')

sdi <- fread('<<<< FILEPATH REDACTED >>>>')[,c('location_id','year_id','sdi'),with=FALSE]
setnames(sdi,c('location_id','year_id'),c('loc_id','year_entering'))
sdi[, year_entering := year_entering - base_year]

message(' ... SDI information sucessfully loaded')

##############################################################################################


##############################################################################################
message(' \nDATA SUBSET TO REGION ONLY ....')

tic('data_prep')
message('Subsetting to region .... \n')
d <- subset(data, region == reg)

# clean up data
message('Cleaning up data .... \n')
d <- clean_dhs(d)
d <- subset(d,childage>=0)
d$yrborn_orig <- d$yrborn

##############################################################################################


##############################################################################################
message('\nCOVARIATES TO UNIT SCALE ....')

# parse out all the fixed effects variables from the gam formula
vars <- get_fe_from_gamformula(formula)

# centrescale all covariates and save for later
# skips SDI/birthorder here since its not loaded in yet and they dont have wild values
cs   <- get_cs_table(d,vars, skipvars = c('sdi','birthorder')) 
for(v in cs$name) {
  d[[paste0(v,'_orig')]] <- d[[v]]
  if(v %in% colnames(d)) 
    d[[v]] <- centscale(d[[v]],v,reverse=FALSE)
}

##############################################################################################




##############################################################################################
message('\nRESHAPE LONG AND FORMAT FOR DTSA ALONG SPECIFIED AGE BINS .... ')

# use survival package handy tools for this
dl <- data.table(survSplit(Surv(childage, died) ~.,cut=ab_times$tstart,data=d))
dl <- merge(dl,ab_times,by='tstart') # add on age bin information
dl[, year_entering := round(yrborn_orig + tstart/12,0)] # identify the year entering each bin
dl <- merge(dl,sdi,by=c('loc_id','year_entering'),all.x=TRUE) # add SDI
dl[,dum := 1] # make a dummy which we use later in period aggregation steps
dl[,fyrint := floor(yrint),by=nid] # floor the year of interview

# survSplit keeps in the last bin which in some cases needs to be censored if the person has not,
#  or would not have (in the case of death) completed the bin yet. Keeping these would create a
#  down bias in the data in recent years due to inflated denominators. As such they are dropped here. 
dl[, cens := (yrborn_orig + tstart/12) + (tstop - tstart)/12 > (cmc_as_year(interview_date_cmc) - base_year) ]
message(sprintf('%i are censored (%f per cent)',sum(dl$cens,na.rm=TRUE),round(sum(dl$cens,na.rm=TRUE)/nrow(dl)*100,2)))
message('Proportion censored by survey nid:')
print(dl[,.(prop_to_censor = sum(cens)/.N), by = nid])
dl <- subset(dl, cens == FALSE)


# anything used as a random effect should be labeled as a factor
dl$mother_id    <- as.factor(dl$mother_id)
dl[, country_ab := as.factor(paste0(country,'_',ab))]
dl$country_orig <- dl$country
dl$country      <- as.factor(dl$country)
dl$nid_yrborn   <- as.factor(paste0(dl$nid,'_',round(dl$yrborn,2)))


# Backfill pre 1970
sdib <- dl[,.(sdib = min(sdi,na.rm=TRUE)),by=loc_id]
dl   <- merge(dl,sdib,by='loc_id')
dl$sdi[is.na(dl$sdi)] <- dl$sdib[is.na(dl$sdi)]


# scale weights so they sum to the same sample size, because of the way mgcv uses weights
# basically the total N should still be the same so we dont change the total contribution to the likelihood
dl[, weight := weight/sum(weight,na.rm=TRUE)*(.N), by = nid] 


toc(log = TRUE) # end data_prep tic

##############################################################################################


##############################################################################################
message(sprintf('\nRUNNING GAM (mgcv::bam) MODEL ON %i CORES AND SAVING OBJECT... ',cores))

# build cluster with X cores and run a bam model
tic('model_fit')
cl  <- makeCluster(cores)
mod <- mgcv::bam(as.formula(formula), data = dl, family = 'binomial', cluster = cl, weight = weight)
stopCluster(cl)

toc(log = TRUE) # end model_fit tic

# print out model summary for the log
print(summary(mod)) 

##############################################################################################


##############################################################################################
message('SAVING FILES .... \n')

# model object and cs description
saveRDS(mod,  file = sprintf('<<<< FILEPATH REDACTED >>>>',rdc_dir))
write.csv(cs, file = sprintf('<<<< FILEPATH REDACTED >>>>',rdc_dir))

# timing logs
tl <- unlist(tic.log(format=TRUE))
tl <- t(matrix(unlist(strsplit(tl,': ')),nrow=2))
write.csv(tl, file = sprintf('<<<< FILEPATH REDACTED >>>>',rdc_dir))

message('SAYING GOODBYE .... \n')
##############################################################################################

