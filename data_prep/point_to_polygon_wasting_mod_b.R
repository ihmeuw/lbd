############### SETUP ###########################

# initiate timing
start_time <- Sys.time()

root <- ifelse(Sys.info()[1]=="Windows", "<<<< FILEPATH REDACTED >>>>>", "<<<< FILEPATH REDACTED >>>>>")
## Load libraries and miscellaneous MBG project functions.

repo <- '<<<< FILEPATH REDACTED >>>>>'
setwd(repo)

# set directories
data_dir <- <<<< FILEPATH REDACTED >>>>>
mbg_dir  <- <<<< FILEPATH REDACTED >>>>>
save_dir <- <<<< FILEPATH REDACTED >>>>>
root_dir <- <<<< FILEPATH REDACTED >>>>>

# Loading libraries, scripts, etc.
package_lib <- ifelse(grepl("geos", Sys.info()[4]),
                      paste0(root,'<<<< FILEPATH REDACTED >>>>>'),
                      paste0(root,'<<<< FILEPATH REDACTED >>>>>'))
.libPaths(package_lib)                                  

source('mbg_central/mbg_functions.R')                   
source('mbg_central/prep_functions.R')                  
source('mbg_central/covariate_functions.R')             
source('mbg_central/misc_functions.R')                  
source('mbg_central/post_estimation_functions.R')
source('mbg_central/gbd_functions.R')
source('mbg_central/shiny_functions.R')
source('mbg_central/holdout_functions.R')
source('mbg_central/polygon_functions.R')
source('mbg_central/collapse_functions.R')
source('mbg_central/graph_data_coverage.R')
source('mbg_central/seegMBG_transform_functions.R')     

package_list <- c('survey', 'pbapply', 'readstata13', 'foreign',
                  'rgeos', 'data.table','raster','rgdal','INLA',
                  'seegSDM','seegMBG','plyr','dplyr', 'foreach',
                  'doParallel', 'bit64', 'mgcv', 'data.table')
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}

## Load iso3 to country name map
iso3_to_country <- fread(paste0(root_dir,"<<<< FILEPATH REDACTED >>>>>"))

indicator <- "wasting_mod_b"

##### 1. ################################################################################################################
## Read all data - this looks like output from ubCov that has been mapped to geographies 
########################################################################################################
# DHS/MICS
all_data <- read.csv(paste0(data_dir,"input_data.csv"), stringsAsFactors = FALSE)
all_data <- as.data.table(all_data)


## fix female coding
all_data$sex[which(all_data$sex == 2)] <- 0

# Somalia data from Damaris Kinyoki
somalia <- read.csv(paste0(root,"<<<< FILEPATH REDACTED >>>>>"))
names(somalia)[names(somalia)=="Latitude"] <- "latitude"
names(somalia)[names(somalia)=="Longitude"] <- "longitude"
names(somalia)[names(somalia)=="Year.of.survey"] <- "start_year"
names(somalia)[names(somalia)=="Number.of.children.examined"] <- "N"
names(somalia)[names(somalia)=="Number.wasting"] <- "wasting_mod_b"
somalia$end_year <- somalia$start_year
somalia$source <- "FSNAU"
somalia$country <- "SOM"
somalia$psu <- 1:nrow(somalia)
somalia$pweight <- NA
somalia$NID <- 270669
somalia$Proportion.wasting <- NULL
somalia$Proportion.stunting <- NULL
somalia$Number.stunting <- NULL
somalia$Proportion.underweight <- NULL
somalia$Number.underweight <- NULL
somalia$Source_Citation <- NULL


##### 2. ################################################################################################################
## Save for data coverage plot
coverage_data <- all_data[, indicator := 1]
coverage_data <- coverage_data[, N := 1]
coverage_data <- coverage_data[, list(N = sum(N)),
                               by = c('indicator','start_year',
                                      'source','latitude','longitude',
                                      'location_code','shapefile',
                                      'country','psu')]


##### 3. ################################################################################################################
## Process your outcome. The result should still be the microdata (rows are individuals)
##    but with a single column for your outcome, i.e. "stunting_binomial" = 0 or 1.
## Example: calculate z-score and whether each individual is stunted (stunting_binomial==1) or not (stunting_binomial==0)

## bring in growth charts
HAZ_chart_months <- read.csv(paste0(mbg_dir,"growth_standards/HAZ_0_60_months.csv"), header = TRUE, sep =",")
WHZ_chart_months <- read.csv(paste0(mbg_dir,"growth_standards/WHZ_0_60_months.csv"), header = TRUE, sep =",")
WAZ_chart_months <- read.csv(paste0(mbg_dir,"growth_standards/WAZ_0_60_months.csv"), header = TRUE, sep =",")

HAZ_chart_weeks <- read.csv(paste0(mbg_dir,"growth_standards/HAZ_0_13_weeks.csv"), header = TRUE, sep =",")
WAZ_chart_weeks <- read.csv(paste0(mbg_dir,"growth_standards/WAZ_0_13_weeks.csv"), header = TRUE, sep =",")

######################################################################################################################
## WHZ

all_data_WHZ <- all_data

# prep all_data -- this needs to happen within each section
all_data_WHZ <- subset(all_data_WHZ, !is.na(age_wks))
all_data_WHZ <- subset(all_data_WHZ, !is.na(child_weight))
all_data_WHZ <- subset(all_data_WHZ, !is.na(child_height))
all_data_WHZ <- subset(all_data_WHZ, all_data_WHZ$child_weight < 99)
all_data_WHZ <- subset(all_data_WHZ, all_data_WHZ$child_height < 999)

# prep all_data_WHZ
all_data_WHZ$child_height <- round_any(all_data_WHZ$child_height, .5)

#remove children whose height is under 45 cm or over 120 cm; the growth charts don't have those values
all_data_WHZ <- all_data_WHZ[all_data_WHZ$child_height >= 45]
all_data_WHZ <- all_data_WHZ[all_data_WHZ$child_height <= 120]

# prep WHZ chart to be joined on
names(WHZ_chart_months)[names(WHZ_chart_months)=="l"] <- "WHZ_l"
names(WHZ_chart_months)[names(WHZ_chart_months)=="m"] <- "WHZ_m"
names(WHZ_chart_months)[names(WHZ_chart_months)=="s"] <- "WHZ_s"
names(WHZ_chart_months)[names(WHZ_chart_months)=="age_cat"] <- "age_cat_1"
names(WHZ_chart_months)[names(WHZ_chart_months)=="length"] <- "child_height"

# merge on growth chart
all_data_WHZ <- merge(all_data_WHZ, WHZ_chart_months, by=c("age_cat_1", "sex", "child_height"), all.x = TRUE, allow.cartesian = TRUE)

# calculate WHZ score
all_data_WHZ$WHZ <- (((all_data_WHZ$child_weight/all_data_WHZ$WHZ_m) ^ all_data_WHZ$WHZ_l)-1)/(all_data_WHZ$WHZ_s*all_data_WHZ$WHZ_l)

## drop unacceptable z scores # https://peerj.com/articles/380/ Crowe, Seal, Grijalva-Eternod, Kerac 2014
all_data_WHZ <- subset(all_data_WHZ, all_data_WHZ$WHZ > -5)
all_data_WHZ <- subset(all_data_WHZ, all_data_WHZ$WHZ < 5)

#################################
## seasonality bias correction ##
#################################

## first, fit the gamms by region
## to do this I need to split dataset into regions
table_file <- paste0(root, "<<<< FILEPATH REDACTED >>>>>")
gaul_table <- read.csv(table_file) %>% data.table

## merge onto df
dat <- base::merge(all_data_WHZ, gaul_table, by.x = "country", by.y = "iso3", all.x = TRUE)
dat$region <- NA

## grab region gaul codes
cssa  <- c(8,49,59,68,76,89)
essa  <- c(43,58,70,77,79,133,150,152,170,205,226,74,257,253,270,40764)
name  <- c(4,40762,40765,145,169,6,248)
sssa  <- c(35,142,172,227,235,271)
wssa  <- c(29,42,45,47,50,66,90,94,106,105,144,155,159,181,182,214,217,221,243)

dat$region[which(dat$gaul %in% cssa)] <- "cssa"
dat$region[which(dat$gaul %in% essa)] <- "essa"
dat$region[which(dat$gaul %in% name)] <- "name"
dat$region[which(dat$gaul %in% sssa)] <- "sssa"
dat$region[which(dat$gaul %in% wssa)] <- "wssa"

## make a time since minimum date
dat$int_year[which(dat$int_year == 93)] <- 1993 ## fix some '93 entries
dat$mo_ind <- dat$int_month + 12 * (dat$int_year - min(dat$int_year, na.rm = TRUE))

## now we can loop through the regions and fit the different seasonal trend
regions <- c("essa", "wssa", "sssa", "name", "cssa")
require(mgcv)
for(reg in regions){
  message(sprintf("On region: %s", reg))
  reg.dat <- subset(dat, region == reg)

  ## fit with uncorrelated errors
  m0 <- gamm(WHZ ~ s(int_month, bs = "cc", k = 12) + s(mo_ind, k = 4) + as.factor(country),
             data = reg.dat)

  ## plot the two components 
  per.plot = TRUE
  if(per.plot){
    png(sprintf("<<<< FILEPATH REDACTED >>>>>/season_gamms_%s_k_4.png", reg),
        width = 20, height = 8, units = "in", res = 300)
    layout(matrix(1:2, ncol = 2))
    plot(m0$gam, scale = 0, main = reg, xlab = "interview month", ylab = "WHZ seasonal peiodic bias")
    dev.off()
    layout(1)
  }

  ## save the model fit
  assign(sprintf("m0_%s", reg), m0)
}

## next, make WHZ adjustments for seasons
## to do this we predict our gamm using one year and one day

reg_periods <- matrix(ncol = 6, nrow = length(seq(1, 12.9, by = .1)))
colnames(reg_periods) <- c("month", regions)
reg_periods[, 1] <- seq(1, 12.9, by = .1)

for(reg in regions){

  ## get a country in the region
  reg.dat <- subset(dat, region == reg)
  ct <- reg.dat$country[!is.na(reg.dat$country)][1]

  pdat <- data.frame(mo_ind    = 1,
                     country   = as.factor(ct), ## random country
                     int_month = seq(1, 12.9, by = .1)) ## random point in time

  reg_periods[, (which(regions %in% reg) + 1)] <- predict(get(sprintf("m0_%s", reg))$gam, newdata = pdat, type = "terms", se.fit = TRUE)$fit[, 2]
}

## now we adjust by each datapoint by the month it was in
dat$WHZ_seas <- NA
reg_periods <- as.data.frame(reg_periods)
total.r <- 0
for(reg in regions){
  for(mo in 1:12){
    rows.to.adjust <- which(dat$region == reg & dat$int_month == mo & !is.na(dat$WHZ))

    total.r <- total.r + length(rows.to.adjust)
    ## find the adjustment
    mean.period <- mean(reg_periods[[reg]]) ## we adjust to the day(s) that's at the mean of the period
    month.val <- reg_periods[[reg]][which(reg_periods[, 1] == mo)]
    delta.adjust <- month.val - mean.period

    ## do the adjustment
    dat$WHZ_seas[rows.to.adjust] <- dat$WHZ[rows.to.adjust] - delta.adjust

  }
}

## check out the differences
adj.r <- which(!is.na(dat$WHZ_seas))
dat$WHZ_seas_complete <- dat$WHZ
dat$WHZ_seas_complete[adj.r] <- dat$WHZ_seas[adj.r]
dat$prev_dif <- (dat$WHZ_seas < -2) - (dat$WHZ < -2) ## difference in prevalence space
dif_seas    <- aggregate(WHZ_seas ~ country + start_year, dat, mean)
dif_no_seas <- aggregate(WHZ      ~ country + start_year, dat, mean)

dif.table <- data.frame(CTRY = dif_no_seas[, 1],
                        year = dif_seas[, 2], 
                        seas = dif_seas[, 3],
                        no_seas = dif_no_seas[, 3],
                        dif = (dif_seas[, 3] - dif_no_seas[, 3]),
                        prev_dif = (dif_seas[, 3] - dif_no_seas[, 3]))

## set the season adjusted data to be WHZ
dat$WHZ <- dat$WHZ_seas
all_data_WHZ <- dat


##################################################################
#### end seasonality correction
##################################################################


## create binary for wasting
all_data_WHZ$wasting_mod_b <- ifelse(all_data_WHZ$WHZ <= -2, 1, 0)
all_data_WHZ$N <- 1

# drop if WHZ is blank
all_data_WHZ <- subset(all_data_WHZ, !is.na(all_data_WHZ$WHZ))

all_data <- as.data.table(all_data_WHZ)
all_data <- all_data[, indicator := 1]

##### 5. ################################################################################################################
## Split up into point and polygon datasets
point_data <- all_data[point==1, ]
poly_data <- all_data[point==0, ]

## Add Somalia dataset
somalia <- as.data.table(somalia)
point_data <- rbind(point_data, somalia, fill = TRUE)

##### 4. ################################################################################################################
## Save for data coverage plot
coverage_data <- rbind(point_data, poly_data, fill = TRUE)
coverage_data <- coverage_data[, list(N = sum(N)),
                               by = c('indicator','start_year',
                                      'source','latitude','longitude',
                                      'location_code','shapefile',
                                      'country','psu')]

##### 5. ################################################################################################################
## Process point_data as you normally would, collapsing to cluster means. Let's call this new dt point_data_collapsed.
all_point_data <- point_data[, list(N=sum(N), wasting_mod_b=sum(wasting_mod_b)),
                             by=c('source', 'start_year','latitude','longitude','country','nid')]
all_point_data$point <- 1

##### 6. ################################################################################################################
## Process poly_data
setnames(all_point_data, "source", "survey_series")
setnames(poly_data, "source", "survey_series")

poly_data_test <- copy(poly_data)
poly_data_test[, N := 1]
poly_data_test <- poly_data_test[, list(N=sum(N)), by=c('start_year', 'country', 'location_code', 'shapefile', 'survey_series')]
poly_data_bad <- poly_data_test[N==1, ]
if(length(poly_data_bad[, survey_series]) > 0) {
    message("This many polygons have 1 observation so will be dropped:")
    print(table(poly_data_bad[, survey_series], poly_data_bad[, start_year]))
##    poly_data_test[, N:= NULL]
    poly_data <- merge(poly_data, poly_data_test, by=c('start_year', 'country', 'location_code', 'shapefile', 'survey_series'))
    poly_data <- poly_data[N.x != N.y, ] ## n.x and n.y are equal where both are 1, i.e. where poly had one cluster
    setnames(poly_data, 'N.x', 'N') ## set the original N col back
    poly_data[, N.y := NULL] ## remove the summed N col
}

## drop strata that have missing pweight so we don't need to drop the whole NID
by_vars <- c('start_year', 'country', 'location_code', 'shapefile', 'survey_series', 'nid')
na.strata <- aggregate(is.na(pweight) ~ start_year + country +
                       location_code + shapefile + survey_series + nid,
                       data = poly_data, sum)

if(sum(na.strata[, 'is.na(pweight)']) > 0){ ## need to drop some
  drop.strata <- na.strata[which(na.strata[, 'is.na(pweight)'] > 0), ]
  drop.rows <- NULL
  for(ds in 1:nrow(drop.strata)){
    drop.rows <- c(drop.rows, which(poly_data$start_year == drop.strata$start_year[ds] &
                                    poly_data$country == drop.strata$country[ds] &
                                    poly_data$location_code == drop.strata$location_code[ds] &
                                    poly_data$shapefile == drop.strata$shapefile[ds] &
                                    poly_data$survey_series == drop.strata$survey_series[ds] &
                                    poly_data$nid == drop.strata$nid[ds])
                   )
  }
  poly_data <- poly_data[-drop.rows, ]
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

point.keepers <- point_data[, c('source', 'start_year','latitude',
                                'longitude','country','nid',
                                "geospatial_id", "master.ind",
                                "cluster_number"), with = FALSE]

na.pw <- aggregate(is.na(pweight) ~ nid, data = poly_data, sum)
drop.nids <- na.pw$nid[which(na.pw[, 2] > 0)]
poly.keepers <- subset(poly_data, !(nid %in% drop.nids))
setnames(poly.keepers, "survey_series", "source")
poly.keepers <- poly.keepers[, c('source', 'start_year','latitude',
                                 'longitude','country','nid',
                                 "geospatial_id", "master.ind",
                                 "cluster_number"), with = FALSE]

keeper.dat <- rbind(point.keepers, poly.keepers)

write.csv(keeper.dat, file = paste0("<<<< FILEPATH REDACTED >>>>>/pre_collapse_wasting_mod_b.csv"), row.names=FALSE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

poly_surveys <- unique(poly_data[, nid])

collapse_each_nid <- function(this_nid) {

    message(paste0('Collapsing NID: ', this_nid))
    test_poly_data <- poly_data[nid==this_nid,]
    test_poly_data$strata <- 0
    names(test_poly_data)[names(test_poly_data)=='cluster_number'] <- 'psu'
    names(test_poly_data)[names(test_poly_data)=='weight'] <- 'pweight'

    # Check for missings
    if(length(test_poly_data$pweight[is.na(test_poly_data$pweight)])>0) {
        message(paste0(length(test_poly_data$pweight[is.na(test_poly_data$pweight)]), ' / ', length(test_poly_data$pweight), ' are missing pweight'))
        return(NULL)
    }

    else {
        collapse_polys <- function(x) {
            setup_design(df = test_poly_data, var = x)
            by_vars <- c('start_year', 'country', 'location_code', 'shapefile', 'survey_series', 'nid')
            keep_vars <- c('geospatial_id', 'master.ind', 'cluster_number')
            poly <- collapse_by(df = test_poly_data,
                                var = x,
                                by_vars = by_vars)
            collapsed <- poly[, c(by_vars, 'mean', 'ss')] ## this converts to ratespace!!
            names(collapsed)[names(collapsed)=='mean'] <- x
            names(collapsed)[names(collapsed)=='ss'] <- 'N'
            collapsed[, eval(x)] <- collapsed[, eval(x)] * collapsed[, 'N'] ## convet back to count space!!
            return(collapsed)
            }
        polys <- c('wasting_mod_b')
        polys <- lapply(polys, collapse_polys)
        merged_polys <- Reduce(function(...) merge(..., all=T), polys)
        return(merged_polys)
    }
}

poly_nids <- unique(poly_data[, nid])
poly_data_collapsed  <- lapply(poly_nids, collapse_each_nid)

# Append all collapsed polygon surveys together
poly_data_collapsed <- do.call(rbind.fill, poly_data_collapsed)
poly_data_collapsed$point <- 0
all_poly_data <- poly_data_collapsed

##### 7. ################################################################################################################
## Append and save a copy for data coverage Shiny (and Jon's static function) before resampling polygons
collapsed <- rbind(all_poly_data, all_point_data, fill=TRUE)

## save the collapsed data so we can report how many point and
## polygons we use from which NIDs post-cleaning
write.csv(collapsed, file = paste0("<<<< FILEPATH REDACTED >>>>>/pre_poly_resampling_wasting_mod_b.csv"), row.names=FALSE)


# ##### 7.1 ################################################################################################################
# ## grab the data for the coverage plots - we plot later to speed up the data processing
coverage_data <- copy(collapsed)
coverage_data <- coverage_data[(shapefile!="COL_DHS_2005" & shapefile!="GTM_DHS_1987") | is.na(shapefile), ]
coverage_data <- coverage_data[, latitude := as.numeric(latitude)]
coverage_data <- coverage_data[, longitude := as.numeric(longitude)]
coverage_data <- coverage_data[, wasting_mod_b := wasting_mod_b / N] #transform from count space to rate space
setnames(coverage_data, 'nid', 'svy_id')
setnames(coverage_data, 'survey_series', 'source')


##### 8. ################################################################################################################
## Resample collapsed polygon data to weighted point data
poly_data_collapsed <- as.data.table(poly_data_collapsed)
poly_data_collapsed$wasting_mod_b_count <- poly_data_collapsed$wasting_mod_b
resampled_poly_data <- resample_polygons_dev(data = poly_data_collapsed,
                                             cores = 30,
                                             indic = 'wasting_mod_b_count') #outcome in count space

all_point_data <- point_data[, list(N=sum(N), wasting_mod_b=sum(wasting_mod_b)), by=c('source', 'start_year','latitude','longitude','country', 'nid')]
all_point_data$point <- 1
all_point_data <- all_point_data[, pseudocluster := FALSE]
all_point_data <- all_point_data[, weight := 1]
all_point_data <- all_point_data[, shapefile := ""]
all_point_data <- all_point_data[, location_code := ""]

resampled_poly_data <- resampled_poly_data[, wasting_mod_b := wasting_mod_b_count]
resampled_poly_data <- resampled_poly_data[, wasting_mod_b_count := NULL]

## rename survey_series to source for poly data
setnames(resampled_poly_data, 'survey_series', 'source')

##### 9. ################################################################################################################
## Append point and polygon collapsed data
all_processed_data <- rbind(all_point_data, resampled_poly_data)
setnames(all_processed_data, 'start_year', 'year')
all_collapsed <- all_processed_data

  ## Replace year with period 1998-2002, 2003-2007, 2008-2012, 2013-2017
  all_collapsed <- subset(all_collapsed, year >= 1997)
  names(all_collapsed)[names(all_collapsed) == "year"] = "original_year"
  all_collapsed <- all_collapsed[original_year >= 1998 & original_year <= 2002, year := 2000]
  all_collapsed <- all_collapsed[original_year >= 2003 & original_year <= 2007, year := 2005]
  all_collapsed <- all_collapsed[original_year >= 2008 & original_year <= 2012, year := 2010]
  all_collapsed <- all_collapsed[original_year >= 2013 & original_year <= 2017, year := 2015]

  all_collapsed <- all_collapsed[, latitude := as.numeric(latitude)]
  all_collapsed <- all_collapsed[, longitude := as.numeric(longitude)]
  all_collapsed <- all_collapsed[!is.na(latitude)]
  all_collapsed <- all_collapsed[!is.na(longitude)]
  all_collapsed <- all_collapsed[latitude>=-90 & latitude<=90]
  all_collapsed <- all_collapsed[longitude>=-180 & longitude<=180]
  all_collapsed <- all_collapsed[, wasting_mod_b := round(wasting_mod_b, 0)]
  ## In clusters where LRI > N (due to tiny samples and every child having LRI), cap at N
  all_collapsed <- all_collapsed[wasting_mod_b > N, wasting_mod_b := N]

write.csv(all_collapsed, file = paste0("<<<< FILEPATH REDACTED >>>>>/wasting_mod_b.csv"), row.names = FALSE)



## 10  ############################################################################################
## now we make the plots once the data has been updated

message("starting data_coverage_plot")
print("starting data_coverage_plot")

coverage_maps <- graph_data_coverage_values(df = coverage_data,
                                            var = 'wasting_mod_b',
                                            title = '',
                                            year_min = '1998',
                                            year_max = '2016',
                                            year_var = 'start_year',
                                            region = 'africa',
                                            sum_by = 'n',
                                            cores = 10,
                                            indicator = 'wasting_mod_b',
                                            high_is_bad = TRUE,
                                            return_maps = TRUE,
                                            legend_title = 'Prevalence \n of MSW')

Sys.time() - start_time
