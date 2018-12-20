## Set repository cloned from Github.
  repo <- <<<< FILEPATH REDACTED >>>>
  root <- ifelse(Sys.info()[1]=="Windows", "<<<<< FILEPATH REDACTED >>>>>", "<<<<< FILEPATH REDACTED >>>>>")

## Load libraries and miscellaneous MBG project functions.
  setwd(repo)
  root <- ifelse(Sys.info()[1]=="Windows", "<<<<< FILEPATH REDACTED >>>>>", "<<<<< FILEPATH REDACTED >>>>>")
  package_lib <- ifelse(grepl("geos", Sys.info()[4]),
                        <<<< FILEPATH REDACTED >>>>,
                        <<<< FILEPATH REDACTED >>>>)
  
  .libPaths(package_lib)                                  # Ensures packages look for dependencies here when called with library(). Necessary for seeg libraries.
  source('mbg_central/mbg_functions.R')                   # Functions to run MBG model.
  source('mbg_central/prep_functions.R')                  # Functions to setup MBG run
  source('mbg_central/covariate_functions.R')             # Functions to prep and transform 5*5 covariates
  source('mbg_central/misc_functions.R')                  # Other computational MBG-related functions.
  source('mbg_central/post_estimation_functions.R')
  source('mbg_central/gbd_functions.R')
  source('mbg_central/shiny_functions.R')
  source('mbg_central/holdout_functions.R')
  source('mbg_central/polygon_functions.R')
  source('mbg_central/collapse_functions.R')
  source('mbg_central/seegMBG_transform_functions.R')     
  package_list <- c('SDMTools', 'survey', 'pbapply', 'readstata13', 'foreign', 'rgeos',
                    'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','dplyr', 'foreach', 'doParallel')
  for(package in package_list) {
    library(package, lib.loc = package_lib, character.only=TRUE)
  }

## Process arguments passed from job submission.
  indicator_name <- as.character(commandArgs()[4])
  target_sex <- as.numeric(commandArgs()[5])
  age_min <- as.numeric(commandArgs()[6])
  age_max <- as.numeric(commandArgs()[7])

## Example arguments
  # indicator_name <- 'edu_mean'
  # target_sex <- 2
  # age_min <- 15
  # age_max <- 49

################################### DATA EXTRACTION #############################################
# Load DHS geography codebook.
    geo_file <- <<<< FILEPATH REDACTED >>>>
    geo_map <- fread(geo_file, verbose = FALSE)
    setnames(geo_map, 'geospatial_id', 'cluster_number')
    geo_map <- geo_map[, cluster_number := as.numeric(cluster_number)]
    geo_point_map <- geo_map[geo_map$point==1]
        dhs_point_map <- unique(geo_point_map[, c('iso3','start_year','cluster_number','lat','long'), with = FALSE])
        dhs_point_map <- dhs_point_map[, type := "DHS"]
        setnames(dhs_point_map, 'start_year', 'year')
    geo_poly_map <- geo_map[geo_map$point!=1]
        dhs_poly_map <- unique(geo_poly_map[, c('iso3','start_year','cluster_number','location_code','admin_level', 'shapefile'), with = FALSE])
        dhs_poly_map <- dhs_poly_map[, type := "DHS"]
        setnames(dhs_poly_map, 'start_year', 'year')
# Load MICS geography codebook.
    geo_file <- <<<< FILEPATH REDACTED >>>>
    geo_map <- fread(geo_file, verbose = FALSE)
    setnames(geo_map, 'geospatial_id', 'cluster_number')
    geo_map <- geo_map[, cluster_number := as.numeric(cluster_number)]
    geo_point_map <- geo_map[geo_map$point==1]
        mics_point_map <- unique(geo_point_map[, c('iso3','start_year','cluster_number','lat','long'), with = FALSE])
        mics_point_map <- mics_point_map[, type := "MICS"]   
        setnames(mics_point_map, 'start_year', 'year')
    geo_poly_map <- geo_map[geo_map$point!=1]
        mics_poly_map <- unique(geo_poly_map[, c('iso3','start_year','cluster_number','location_code','admin_level', 'shapefile'), with = FALSE])
        mics_poly_map <- mics_poly_map[, type := "MICS"]
        setnames(mics_poly_map, 'start_year', 'year')
        
## Load all single-year point data for DHS/MICS (data where education attainment is coded continuously).
    data <- fread(<<<< FILEPATH REDACTED >>>>)

## Load iso3 map.
    dhs_mics_iso3s <- list.files(<<<< FILEPATH REDACTED >>>>, pattern = 'codebook_survey_gradeLevel', full.names = TRUE)
    dhs_mics_iso3s <- rbindlist(lapply(dhs_mics_iso3s, fread))
    dhs_mics_iso3s <- unique(dhs_mics_iso3s[, c('nid','iso3'), with = FALSE])
    dhs_iso3s <- fread(<<<< FILEPATH REDACTED >>>>)
    dhs_iso3s <- unique(dhs_iso3s[, c('nid','iso3'), with = FALSE])
    dhs_mics_iso3s <- rbind(dhs_mics_iso3s, dhs_iso3s)
    dhs_mics_iso3s <- dhs_mics_iso3s[, nid := gsub('y','',nid)]
    data <- merge(data, dhs_mics_iso3s, by='nid')

## Get correct subset.
    subset_data <- data[age_start >= age_min & age_start <= age_max & sex == target_sex, ]
    subset_data <- subset_data[, total := NULL]
    subset_data <- subset_data[, total_ss := NULL]

# Mean and standard error of mean for points.
    point_means <- subset_data[,.(total=sum(count),
                                  total_ss=sum(sample_size),
                                  mean=weighted.mean(x=edu_yrs,w=proportion),
                                  sd=wt.sd(x=edu_yrs,w=proportion)),
                               by=.(year,iso3,type,ihme_loc_id,nid)]
    point_means <- point_means[,mean_se:=sd / sqrt(total_ss)]
    setnames(point_means, 'ihme_loc_id', 'cluster_number')
    point_means <- point_means[, cluster_number := as.numeric(cluster_number)]

## Merge maps and only keep mapped NIDs. This is all point data for DHS/MICS.
    # DHS
        dhs_point_means <- merge(point_means, dhs_point_map, by=c('type', 'iso3','year','cluster_number'))
    # MICS
        mics_point_means <- merge(point_means, mics_point_map, by=c('type', 'iso3','year','cluster_number'))
    # Append all points
        all_dhs_mics_point_means <- rbind(dhs_point_means, mics_point_means)

## Merge on location_code/shapefile maps to DHS/MICS points and collapse the ones that merge to polygons.
    poly_subset_data <- copy(subset_data)
    setnames(poly_subset_data, 'ihme_loc_id', 'cluster_number')
    poly_subset_data <- poly_subset_data[, cluster_number := as.numeric(cluster_number)]
    dhs_poly_subset_data <- merge(poly_subset_data, dhs_poly_map, by=c('type','iso3','year','cluster_number'))
    mics_poly_subset_data <- merge(poly_subset_data, mics_poly_map, by=c('type','iso3','year','cluster_number'))
    poly_subset_data <- rbind(dhs_poly_subset_data,
                              mics_poly_subset_data)

# Mean and standard error of mean for polygons.
    poly_means <- poly_subset_data[,.(total=sum(count),
                                      total_ss=sum(sample_size),
                                      mean=weighted.mean(x=edu_yrs,w=proportion),
                                      sd=wt.sd(x=edu_yrs,w=proportion)),
                                   by=.(year,iso3,type,location_code,shapefile,nid)]
    dhs_mics_poly_means <- poly_means[,mean_se:=sd / sqrt(total_ss)]

## Load all single-year polygon data for IPUMS.
    ipums_poly_data <- fread(<<<< FILEPATH REDACTED >>>>)

## Map NID to IPUMS sample number to merge on shapefile names.
    ipums_sample_to_nid <- fread(<<<< FILEPATH REDACTED >>>>)
    ipums_sample_to_nid <- ipums_sample_to_nid[, c('iso3','nid'), with=FALSE]
    ipums_poly_data <- ipums_poly_data[, nid := as.numeric(nid)]
    ipums_poly_data <- merge(ipums_poly_data, ipums_sample_to_nid, by='nid')
    # Merge shapefile names from my map by sample number.
    ipums_sample_to_shapefile <- fread(<<<< FILEPATH REDACTED >>>>)
    ipums_sample_to_shapefile <- unique(ipums_sample_to_shapefile[, c('shapefile','iso3','year'), with=FALSE])
    ipums_poly_data <- merge(ipums_poly_data, ipums_sample_to_shapefile, by=c('iso3','year'))

## Get correct IPUMS (census) subset.
    ipums_poly_subset_data <- ipums_poly_data[age_start >= age_min & age_start <= age_max & sex == target_sex, ]
    ipums_poly_subset_data <- ipums_poly_subset_data[, total := NULL]
    ipums_poly_subset_data <- ipums_poly_subset_data[, total_ss := NULL]
    setnames(ipums_poly_subset_data, 'ihme_loc_id', 'location_code')

# Get IPUMS polygon totals.
    ipums_poly_total <- ipums_poly_subset_data[,.(total=sum(count),
                                                  total_ss=sum(sample_size)),
                                               by=.(year,iso3,location_code,shapefile,type,nid)]
# Mean and standard error of mean for polygons.
    ipums_poly_means <- ipums_poly_subset_data[,.(mean=weighted.mean(x=edu_yrs,w=proportion),
                                                  sd=wt.sd(x=edu_yrs,w=proportion)),
                                               by=.(year,iso3,type,location_code,shapefile,nid)]
    ipums_poly_means <- data.table(merge(ipums_poly_means,ipums_poly_total,by=c("year","iso3","location_code","shapefile","type","nid")))
    ipums_poly_means <- ipums_poly_means[,mean_se:=sd / sqrt(total_ss)]

## Process binned point data for DHS/MICS (attainment data coded to bins, such as 0-6 years, 7-12 years, etc.).
    binned_point_dir <- <<<< FILEPATH REDACTED >>>>
    binned_point_files <- list.files(binned_point_dir)
    map_binned_dhs_mics <- function(x) {
        cluster_data <- fread(paste0(binned_point_dir, '/', x))
        iso3 <- strsplit(x, split="_")
        iso3 <- iso3[[1]][3]
        cluster_data <- cluster_data[, iso3 := iso3]
        ## Get correct subset.
        subset_data <- cluster_data[age_start >= age_min & age_start <= age_max & sex == target_sex, ]
        subset_data <- subset_data[, total := NULL]
        subset_data <- subset_data[, total_ss := NULL]
        # Mean and standard error of mean for points.
        point_means <- subset_data[,.(total=sum(count),
                                      total_ss=sum(sample_size),
                                      mean=weighted.mean(x=edu_yrs,w=proportion),
                                      sd=wt.sd(x=edu_yrs,w=proportion)),
                                   by=.(year,iso3,type,ihme_loc_id,nid)]
        point_means <- point_means[,mean_se:=sd / sqrt(total_ss)]
        setnames(point_means, 'ihme_loc_id', 'cluster_number')
        point_means <- point_means[, cluster_number := as.numeric(cluster_number)]
        return(point_means)
    }
    binned_dhs_mics_point_means <- rbindlist(lapply(binned_point_files, map_binned_dhs_mics))

## Map binned points.
    # DHS
        binned_dhs_point_means <- merge(binned_dhs_mics_point_means, dhs_point_map, by=c('iso3','year','type','cluster_number'))
    # MICS
        binned_mics_point_means <- merge(binned_dhs_mics_point_means, mics_point_map, by=c('iso3','year','type','cluster_number'))
    # Append all points
        binned_dhs_mics_point_means <- rbind(binned_dhs_point_means, binned_mics_point_means)

## Process binned polygon data for DHS/MICS.
    parse_survey_names <- function(x) {
        x <- strsplit(x, split="_")
        x <- paste(x[[1]][1],x[[1]][2],x[[1]][3],sep = "_")
        return(x)
    }
    binned_point_files <- list.files(binned_point_dir, full.names = FALSE)
    binned_surveys <- unique(lapply(binned_point_files, parse_survey_names))
    collapse_binned_polys <- function(survey) {
        # Read all clusters associated with this survey
        iso3 <- strsplit(survey, split="_")
        iso3 <- iso3[[1]][3]
        all_clusters <- list.files(binned_point_dir, pattern=survey, full.names = TRUE)
        all_clusters <- rbindlist(lapply(all_clusters, fread))
        setnames(all_clusters, 'ihme_loc_id', 'cluster_number')
        all_clusters <- all_clusters[, cluster_number := as.numeric(cluster_number)]
        all_clusters <- all_clusters[, iso3 := iso3]
        # Merge maps
        if(grepl("DHS", survey)) all_clusters <- merge(all_clusters, dhs_poly_map, by=c('iso3','year','type','cluster_number'))
        if(grepl("MICS", survey)) all_clusters <- merge(all_clusters, mics_poly_map, by=c('iso3','year','type','cluster_number'))
        # Mean and standard error of mean for polygons
        poly_means <- all_clusters[,.(total=sum(count),
                                      total_ss=sum(sample_size),
                                      mean=weighted.mean(x=edu_yrs,w=proportion),
                                      sd=wt.sd(x=edu_yrs,w=proportion)),
                                   by=.(year,iso3,type,location_code,shapefile,nid)]
        poly_means <- poly_means[,mean_se:=sd / sqrt(total_ss)]
        if(length(poly_means$iso3)==0) print(paste0('Survey: ', survey))
        return(poly_means)
            
    }
    binned_dhs_mics_poly_means <- rbindlist(lapply(binned_surveys, collapse_binned_polys))

## Process binned polygon data for IPUMS
    binned_ipums_dir <- <<<< FILEPATH REDACTED >>>>
    binned_poly_files <- list.files(binned_ipums_dir, full.names = FALSE)
    binned_ipums_surveys <- unique(lapply(binned_poly_files, parse_survey_names))
    collapse_binned_ipums_polys <- function(survey) {
        # Read all clusters associated with this survey
        iso3 <- strsplit(survey, split="_")
        iso3 <- iso3[[1]][3]
        all_polys <- list.files(binned_ipums_dir, pattern=survey, full.names = TRUE)
        all_polys <- rbindlist(lapply(all_polys, fread))
        setnames(all_polys, 'ihme_loc_id', 'location_code')
        all_polys <- all_polys[, location_code := as.numeric(location_code)]
        all_polys <- all_polys[, nid := as.numeric(nid)]
        all_polys <- all_polys[, iso3 := iso3]
        # Merge shapefile names from my map by iso3, year, location_code
            ipums_iso3_to_shapefile <- unique(ipums_sample_to_shapefile[, c('shapefile','year','iso3'), with=FALSE])
            ipums_poly_data <- merge(all_polys, ipums_iso3_to_shapefile, by=c("iso3","year"))
        # Subset
            ipums_poly_subset_data <- ipums_poly_data[age_start >= age_min & age_start <= age_max & sex == target_sex, ]
            ipums_poly_subset_data <- ipums_poly_subset_data[, total := NULL]
            ipums_poly_subset_data <- ipums_poly_subset_data[, total_ss := NULL]
        # Mean and standard error of mean for polygons
            ipums_poly_means <- ipums_poly_subset_data[,.(total=sum(count),
                                                          total_ss=sum(sample_size),
                                                          mean=weighted.mean(x=edu_yrs,w=proportion),
                                                          sd=wt.sd(x=edu_yrs,w=proportion)),
                                                       by=.(iso3,year,nid,type,location_code,shapefile)]
            ipums_poly_means <- ipums_poly_means[,mean_se:=sd / sqrt(total_ss)]
        return(ipums_poly_means)
            
    }
    binned_ipums_poly_means <- rbindlist(lapply(binned_ipums_surveys, collapse_binned_ipums_polys))
    
## Append everything. Merge iso3 to each dataset first, because IPUMS binned polys will have iso3.
    dhs_mics_poly_means <- dhs_mics_poly_means[, point := 0]
    dhs_mics_poly_means <- dhs_mics_poly_means[, binned := 0]
    
    binned_dhs_mics_poly_means <- binned_dhs_mics_poly_means[, point := 0]
    binned_dhs_mics_poly_means <- binned_dhs_mics_poly_means[, binned := 1]
    
    all_dhs_mics_point_means <- all_dhs_mics_point_means[, point := 1]
    all_dhs_mics_point_means <- all_dhs_mics_point_means[, binned := 0]
    
    binned_dhs_mics_point_means <- binned_dhs_mics_point_means[, point := 1]
    binned_dhs_mics_point_means <- binned_dhs_mics_point_means[, binned := 1]
    
    ipums_poly_means <- ipums_poly_means[, point := 0]
    ipums_poly_means <- ipums_poly_means[, binned := 0]
    
    binned_ipums_poly_means <- binned_ipums_poly_means[, point := 0]
    binned_ipums_poly_means <- binned_ipums_poly_means[, binned := 1]
    
    all_edu_data <- rbind(dhs_mics_poly_means,
                          binned_dhs_mics_poly_means,
                          all_dhs_mics_point_means,
                          binned_dhs_mics_point_means,
                          ipums_poly_means,
                          binned_ipums_poly_means,
                          fill = TRUE)

    all_edu_data <- all_edu_data[, lat := as.numeric(lat)]
    all_edu_data <- all_edu_data[, long := as.numeric(long)]

## Save for data coverage plots.
    coverage_data <- copy(all_edu_data)
    setnames(coverage_data, 'total_ss', 'N')
    setnames(coverage_data, 'type', 'source')
    setnames(coverage_data, 'cluster_number', 'psu')
    setnames(coverage_data, 'iso3', 'country')
    setnames(coverage_data, 'lat', 'latitude')
    setnames(coverage_data, 'long', 'longitude')
    setnames(coverage_data, 'nid', 'svy_id')
    coverage_data <- coverage_data[, svy_id := gsub('.DTA','',svy_id)]
    coverage_data <- coverage_data[, svy_id := as.numeric(svy_id)]
    coverage_data <- coverage_data[, survey := paste0(source,'_',country,'_',year)]
    ## Test for duplicated polygon ids.
    test_polys <- coverage_data[is.na(latitude),]
    test_polys <- test_polys[!is.na(location_code),]
    for(this_survey in unique(test_polys[, survey])) {
        if(length(unique(test_polys[survey==this_survey, c('location_code','shapefile'), with=F]$location_code)) !=
            length((test_polys[survey==this_survey, c('location_code','shapefile'), with=F]$location_code))) {
            message(this_survey)
        }
    }

## Data coverage plot
  coverage_maps <- graph_data_coverage_values(df = coverage_data,
                                              var = 'mean',
                                              title = 'Education',
                                              year_min = '1998',
                                              year_max = '2016',
                                              year_var = 'year',
                                              region = 'africa',
                                              sum_by = 'n',
                                              cores = 10,
                                              indicator = indicator_name,
                                              high_is_bad = FALSE,
                                              return_maps = TRUE,
                                              legend_title = "Mean years",
                                              save_on_share = TRUE,
                                              simplify_polys = TRUE)

## Resample polygon to points.
  resampled_poly_data <- resample_polygons_dev(data = all_edu_data,
                                               cores = 20,
                                               indic = 'mean')

## Final data formatting.
  final_data <- all_edu_data[point==1, ]
  final_data <- final_data[, pseudocluster := FALSE]
  final_data <- final_data[, weight := 1]
  setnames(final_data, c('lat','long'), c('latitude','longitude'))
  final_data <- rbind(final_data, resampled_poly_data)
  final_data <- subset(final_data, year >= 1998)
  final_data <- final_data[, original_year := year]
  setnames(final_data, 'mean', indicator_name)
  setnames(final_data, 'iso3', 'country')
  setnames(final_data, 'type', 'source')
  setnames(final_data, 'total_ss', 'N')
  
write.csv(final_data, file = <<<< FILEPATH REDACTED >>>>, row.names = FALSE)
