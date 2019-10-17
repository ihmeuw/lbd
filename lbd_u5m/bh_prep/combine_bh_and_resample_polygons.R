## #######################################################
## PURPOSE: Combine prepared SBH and CBH datasets, polygon resample
## INPUTS: in <<<< FILEPATH REDACTED >>>>:
##            prepared_sbh.RDS and prepared_cbh.RDS
## OUTPUTS: one file, resampled to points, for input data
##         This file is the input for mbg modeling
## #######################################################

## #############################################################
## SETUP

# user name
user       <- Sys.info()['user']

# libraries
library(data.table)
library(gtools)
library(parallel)
library(raster)
library(seegMBG)
library(magrittr)

# dirs
mbgrepo <- sprintf('<<<< FILEPATH REDACTED >>>>',user)
ig_repo         <- sprintf('<<<< FILEPATH REDACTED >>>>',user) 
datadir <- '<<<< FILEPATH REDACTED >>>>'

## Load libraries and MBG project functions.
commondir <- sprintf('<<<< FILEPATH REDACTED >>>>')
package_list<-c(t(read.csv(sprintf('<<<< FILEPATH REDACTED >>>>',commondir),header=FALSE)))
source(paste0(mbgrepo, '<<<< FILEPATH REDACTED >>>>'))
mbg_setup(package_list = package_list, repos = mbgrepo)
load_mbg_functions(ig_repo)
## #############################################################

## ##############################################################
## Load in data and merge to one file
## These will have already been prepared and come directly out of the 
## following scripts: cbh_collapse.R and sbh_collapse.R

# cbh
cbh <- readRDS("<<<< FILEPATH REDACTED >>>>")

# sbh
sbh <- readRDS("<<<< FILEPATH REDACTED >>>>")

# combine
bh <- setDT(gtools::smartbind(sbh, cbh))

# clean up some naming things
setnames(bh, c('latnum','longnum'), c('latitude','longitude'))

shapefile_version <- "current"
modeling_shapefile_version <- NULL
## #############################################################


## ##############################################################
## Polygon resample

# pass it a unique geography-year list and simply merge it later on all age bins
#  we will apply the resample to each age bin later
bh_point <- subset(bh, point == 1)
bh_poly  <- subset(bh, point == 0)

tmp <- unique(bh_poly[,c('geo_unit','year','location_code','shapefile','point','latitude','longitude'),with=FALSE])
tmp[, died:=1]; tmp[,N:=1]

#countries we have SBH polygon data for that needs to go through resample
extra_countries <- c("TTO", "ISR", "ALB", "ARM", "AZE", "CHL", "CUB", "GEO", "KAZ", "MDV", "TUR", "UKR")
gaul_codes <-  unique(c(get_adm0_codes(c('stage1','stage2',extra_countries),
                                       shapefile_version = shapefile_version)))

# run the resample
rp <- resample_polygons(
  tmp,
  cores            = 20,
  ignore_warnings  = TRUE,
  indic            = 'died',
  pull_poly_method = "fast",
  density          = 0.0001,
  gaul_list        = gaul_codes,
  shapefile_version = shapefile_version
)

# merge back the resampled points onto the dataset
# the following variables are new: weight, pseudocluster // lat and long are updated
bh_poly[, latitude := NULL]
bh_poly[, longitude := NULL]
bh_poly[,weight := NULL]
bh_point[,pseudocluster := FALSE]
bh_point[,weight := 1]
rp[,died := NULL]
rp[,N := NULL]

# Many to many merge, since rp has many points per geo_unit and bh_poly has several age bins
# We can deal with this using `cartesian=TRUE` in the merge command
bh_rpoly <- merge(
  x   = bh_poly,
  y   = rp,
  by  = c('geo_unit','year','location_code','shapefile','point'),
  allow.cartesian = TRUE
)

# bring back together poly and point for one final point file
bh <- rbind(bh_rpoly, bh_point, fill = T) 


# add age bins labels make sure this matches sbh prep and cbh prep
tmp <- data.table(age=1:length(unique(bh$ab)),ab=unique(bh$ab)) # make sure these are ordered correctly
bh <- merge(bh, tmp, by='ab', all.x=TRUE)
## #############################################################

## ##############################################################
## Use a weighted sample size appraoch rather than a directly 
##  weighed log-likelihood approach as we used before
bh[point == 1, weight := 1]
bh[is.na(sbh_weight), sbh_weight:=1] # sbh weight does not apply for CBH observations
bh[, weight := sbh_weight*weight]

bh[,died:=died*weight]
bh[,N:=N*weight]
bh[,weight:=1]

## ##############################################################

# make a dummy for sbh data, gets used in modelling later
bh[,sbh := as.numeric(data_type == 'sbh')]

# Remove rows with na in died or N
if(sum((is.na(bh$died) | is.na(bh$N))) > 0)
  message(sprintf('WARNING: There are %s rows (%f%%) where died or N is NA.\n',
                  prettyNum(sum((is.na(bh$died) | is.na(bh$N))),big.mark=','),
                  sum((is.na(bh$died) | is.na(bh$N)))/nrow(bh)*100))
bh <- na.omit(bh, cols=c('died','N'))

# archived version
rd <- gsub(" ","_",gsub(":","_",gsub("-","_",Sys.time())))
saveRDS(bh, sprintf('<<<< FILEPATH REDACTED >>>>',datadir,rd))

# ready to model version
saveRDS(bh,sprintf('<<<< FILEPATH REDACTED >>>>',datadir))
