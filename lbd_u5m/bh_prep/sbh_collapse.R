## #######################################################
## PURPOSE: This script combines adjusted SBH yearly data for each
#           SBH survey into one file and collapses onto the desired
#           years run_sbh_prep_models_predict.R.
#  OUPUT: A dated prepped sbh combined and collapsed dataset, and
#         an undated (most up to date one save in <<<< FILEPATH REDACTED >>>>)
## #######################################################




## ##############################################################
## USER DEFINED GLOBALS FOR RUN DATE AND PERIOD TABULATION
user       <- Sys.info()['user']


# Define the sbh prediction run date to source individual prepped surveys from.
sbh_prep_rd <- '<<<< DATE REDACTED >>>>'

# Start of each period. 

period_start <- c(1998:2017)
period_label <- c(1998:2017)

# NOTE: Age bins have already been defined in the model fitting process, 
#       these too much match what is being done in cbh_collapse.R

## ##############################################################






## ##############################################################
## SETUP - LOOK FOR COMPLETED NID PREPS TO BRING IN, WARN ABOUT MISSING ONES

# libraries
library(data.table)
library(parallel)

# root dir for outputs
root <- '<<<< FILEPATH REDACTED >>>>'

# load custom functions
setwd(sprintf('<<<< FILEPATH REDACTED >>>>', user))
source('<<<< FILEPATH REDACTED >>>>')

# load out of sample variance by age bin from the DHS validation work
oos <- fread('<<<< FILEPATH REDACTED >>>>')
oos <- oos[,c('ab','var_logit_lss','var_lss','var_logit_lss_w','var_lss_w','SD_res_EB','cv_logit_lss','cv_res_EB'),with=FALSE] 
oos[, ab := toupper(ab)]
oos[, ab := gsub('Y','yr',ab)]

# source a list of all the NIDs that should be in there:
nid_dirs <-  list.dirs(sprintf('<<<< FILEPATH REDACTED >>>>',root,sbh_prep_rd))[-1]
all_nids <- gsub(sprintf('<<<< FILEPATH REDACTED >>>>',root,sbh_prep_rd),'',nid_dirs)


# look through each directory and search for a completed file, record where we have and dont have data
missing_prep_nids <- missing_prep_dirs <- as.character()
for(d in 1:length(nid_dirs)){
  if(!file.exists(sprintf('%s/fin.txt',nid_dirs[d]))) {
    missing_prep_nids <- c(missing_prep_nids,all_nids[d])
    missing_prep_dirs <- c(missing_prep_dirs,nid_dirs[d])
  }    
}

nids_with_data <- all_nids[!all_nids %in% missing_prep_nids]
dirs_with_data <- nid_dirs[!nid_dirs %in% missing_prep_dirs]

if(length(missing_prep_nids) > 0) {
  message(sprintf('%i EXPECTED SBH PREP FILES ARE MISSING!\n',length(missing_prep_nids)))
  message(paste('THE FOLLOWING NIDS DO NOT HAVE PREP FILES:',paste(missing_prep_nids,collapse=', '),'\n'))
  message('WARNING: Please make sure all predict runs have finished, 
          and if so and there are still missing files, check logs for errors.
          this code will continue to collapse and combine only those surveys 
          which have prepped data for.')
}

## ##############################################################





## ##############################################################
## For those with data, collapse and summarize them

# write a function that loads data, collapses and summarizes
load_collapse_summarize <- function(dir){
  # store the nid
  nid <-  gsub(sprintf('<<<< FILEPATH REDACTED >>>>',root,sbh_prep_rd),'',dir)
  
  # load in the data, each subset file
  files <- list.files(dir, pattern = 'prepped_yearly_aggregation_draws_')
  d <- list()
  for(f in files)
    d[[f]] <- readRDS(sprintf('%s/%s',dir,f))
  d <- do.call('rbind', d)
  
  d[, period := findInterval(year_entering, period_start)]
  
  # keep only those with any estimates in a period greater than 0
  d <- subset(d, period > 0) # could also check based on svyyr
  
  if(nrow(d) == 0){ # if there are no data in the study period, end the function
    
    message(sprintf('FOR NID %s NO OBSERVATIONS WITHIN STUDY PERIOD. RETURNING BLANK DT',nid))
    return(data.table())
    
  } else { # if there are data within the correct period, then continue on with the function
  
    # remove year entering so can c  ollapse on period instead
    d$year_entering <- NULL
    drawcols <- colnames(d)[grep('V\\d*',colnames(d))]
    
    # merge on the oos variance
    d <- merge(d, oos, by = 'ab', all.x = TRUE)

    # Todo: Incorporate uncertainty in EEB as well? 
    
    ## downweight SBH using a simulation approach
    ## NOTE Uncertainty in p is very small, leading to almost the same sample sizes as without this correction. 
    #      For now do not do a correction.. 
    #      Also to note: there are often cases of zero deaths with small sample sizes in this , which  leads to issues
    #      May need to do this at the child level in the prediction step as bernoulli, then we get real SS as well
    # 1. simulate draws of Y using draws of p and EEB. Make extra draws to try and avoid zeroes
    drawmultiplier <- 10
    d2 <- matrix(NA,nrow=nrow(d),ncol=length(drawcols)*drawmultiplier)
  #  d3 <- matrix(NA,nrow=nrow(d),ncol=length(drawcols)*drawmultiplier)
    d4 <- matrix(NA,nrow=nrow(d),ncol=length(drawcols)*drawmultiplier)
    for(dr in 1:length(drawcols))
     for(drset in 0:(drawmultiplier-1)){
       deeb <- ceiling(d$EEB) #ceiling(rnorm(nrow(d2),d$EEB,d$EEB*d$cv_res_EB)) #ceiling(d$EEB) #
       deeb[deeb<0] <- 0 # resolves a few edge cases
       #pppp <- plogis(rnorm(nrow(d),qlogis(d[[drawcols[dr]]]),sqrt(d$var_logit_lss))) # use CV instead?
       pppp <- plogis(rnorm(nrow(d),qlogis(d[[drawcols[dr]]]),qlogis(d[[drawcols[dr]]])*d$cv_logit_lss)) # scaled residual error
       d2[,dr+length(drawcols)*drset] <- rbinom(nrow(d2),size=deeb,prob=pppp)
       d4[,dr+length(drawcols)*drset] <- pppp
       
     #  d3[,dr+length(drawcols)*drset] <- deeb
       #d2[,dr+length(drawcols)*drset] <- rbinom(nrow(d2),size=deeb,prob=d[[drawcols[dr]]]) #/deeb
    }
    ## 2. calculate new EEB as var(Y/N) / mean(Y/N)*(1-mean(Y/N))
    ##    save both old and new EEB
    # small numbers of deaths lead to difficult emprical estimates here. 
    #varofdraws   <- apply(d2,1,var,na.rm=TRUE) 
    #meanofdraws  <- apply(d4,1,mean)  #apply(d2,1,mean,na.rm=TRUE) # # #  #  # apply(d[,drawcols,with=FALSE],1,var) #apply(d2,1,mean) # 
    #newN2      <-  varofdraws/( meanofdraws * (1-meanofdraws) )  #+d$var_lss) # add on out of sample variance
    #s2 <- newN2/d$EEB
    #newN      <- rowMeans(d[,drawcols,with=FALSE]*(1-d[,drawcols,with=FALSE]) / varofdraws)
    #d$sbh_weight <- newN/d$EEB 
    #summary(d$sbh_weight)
    #d[,.(m=mean(sbh_weight)),by=ab]
    
    # do it the Beta Binomial way using method of moments
    n  <- ncol(d2)
    m1 <- apply((d2),1,sum)/n
    m2 <- apply((d2)^2,1,sum)/n
    
    a <- (d$EEB*m1-m2)/( d$EEB * (m2/m1 - m1 - 1) + m1 )
    b  <- ( (d$EEB - m1)*(d$EEB - m2/m1) ) / ( d$EEB * (m2/m1 - m1 -1 ) + m1)
      
    newN  <- ( (a/(a+b)) * (1-(a/(a+b))) * d$EEB * (a+b)^2 * (a+b+1)) / ( a*b*(a+b+d$EEB) )  
    
    #vBB <- (d$EEB*a*b*(a+b+d$EEB)) / ( ( a+b)^2 * (a +b +1) )
    #vB  <- apply(d2,1,var,na.rm=TRUE)
  
    
    s <- newN/d$EEB  
    #summary(s)
    #sum(s>1)
    #data.table(s=s,a=d$ab)[,.(ms=mean(s)),by=a]
    d$sbh_weight <- s 
    d$sbh_logit_variance <- apply(qlogis(d4),1,var,na.rm=TRUE) 
    rm(d4)
    # summary(d$sbh_weight)
     
     # add on some variance components to be used in modelling later, maybe
     #d$sbh_prepmodel_variance       <- apply(d[,drawcols,with=FALSE],1,var,na.rm=TRUE) 
     #d$sbh_prepmodel_logit_variance <- apply(qlogis(as.matrix(d[,drawcols,with=FALSE])),1,var,na.rm=TRUE) 


    # collapse to the geo_unit-period-age bin: weighted mean each draw, sum EEB
    # NOTE, this is expecting draws to be of format VXXXX, V and up to 4 digits.  
    collapse_on  <- drawcols # columns for draw variables, (Must be V then some numbers)
    collapse_by  <- colnames(d)[grep('V\\d*|EEB|sum_wgt',colnames(d),invert=TRUE)]
    d <- d[, c(list(N = sum(EEB), Nwgt = sum(sum_wgt)), lapply(.SD, weighted.mean, sum_wgt)), .SDcols = collapse_on, by = collapse_by]
  
    # summarize draws, mean upper lower, 2.5% and 97.5%
    summ <- data.table(t(d[,apply(.SD,1,quantile,probs=c(0.5,0.025,0.975)),.SDcols=collapse_on]))
    colnames(summ) <- c('haz','haz_lower','haz_upper')
    
    # bind just the summaries back on 
    d <- cbind(d[, c(collapse_by, 'N', 'Nwgt'), with = FALSE], summ)
    
    # return the dataset
    return(d)
  }
}

# apply over each survey and do the thing
message(sprintf('Loading, collapsing, and summarizing %i surveys.', length(dirs_with_data)))
sbh <- mclapply(dirs_with_data, load_collapse_summarize, mc.cores = 20)

# sbh <- NULL
# for(dir in dirs_with_data){
#   d <- load_collapse_summarize(dir)
#   message(sprintf('Done with %s',dir))
#   sbh <- rbind(sbh, d)
# }

# rbind over the list to get all the data into one frame
sbh <- do.call('rbind',sbh)
sbh$year <- period_label[sbh$period]

# Drop rows where year is greater than survey year
# Drop rows where year is greater than survey year
if(sum(sbh$year > sbh$svyyr) > 0)
  message(sprintf('WARNING: There are %s rows (%f%%) where year is greater than survey year.\n',
                  prettyNum(sum(sbh$year > sbh$svyyr),big.mark=','),
                  sum(sbh$year > sbh$svyyr)/nrow(sbh)*100))
sbh <- sbh[year <= svyyr]

# median estimate of died for each row
sbh[, died := N*haz]

# clean up some variables to make sure the final prepped sbh dataset matches the final prepped cbh dataset
sbh[, data_type := 'sbh']
info <- fread('<<<< FILEPATH REDACTED >>>>')[,c('country','loc_id'),with=FALSE]
sbh <- merge(sbh, info, by='loc_id', all.x = TRUE)
sbh[,loc_id := NULL]

# get rid of leading numbers on ab, if they are there.
if(any(grepl('1_', sbh$ab))) sbh[, ab := substring(ab, 3)] 

# print a quick summary of the data
print(sbh[, .(haz = weighted.mean(haz,N), N = sum(N)), by = .(ab,period)][order(ab,period)])
message(sprintf('There are %s approx. births in the SBH dataset', prettyNum(sum(sbh[ab=='NN']$N), big.mark = ',')))
message(sprintf('There are %s approx. deaths in the SBH dataset', prettyNum(sum(sbh$died), big.mark = ',')))
print(summary(sbh$haz_upper))

### ##############################################################


## Cleanup
sbh[, c('var_logit_lss_w', 'var_lss_w', 'SD_res_EB', 'cv_logit_lss', 'cv_res_EB', 'var_logit_lss', 'var_lss', 'Nwgt') := NULL]

# with small N, variance could have been zero so NA sbh_weight. make these one
sbh[is.na(sbh_weight), sbh_weight := 1]
# with small N there are some cases where it is >1 due to MC variability. make it one
sbh[sbh_weight>1, sbh_weight := 1]
sbh[sbh_weight<0, sbh_weight := 1] # super tiny N's fix this though

## ##############################################################
## SAVE

# Save a dated version (with the run_date used for the individual surveys)
saveRDS(sbh, '<<<< FILEPATH REDACTED >>>>')

# Save an undated file as well, the one that other code will reference
saveRDS(sbh, '<<<< FILEPATH REDACTED >>>>')


### ##############################################################





# TODO: Create SBH weights

