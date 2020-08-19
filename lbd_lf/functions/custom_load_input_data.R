## Load input data from required location
#   Arguments:
#     indicator = Specific outcome to be modeled within indicator category, i.e. "edu_0"
#     simple    = Single polygon that defines boundaries of the entire area you want to model over.
#   Returns: Input data subset to modeling area.
load_input_data <- function(indicator, simple = NULL, agebin = 0, removeyemen = FALSE, pathaddin = "",
                            withdate=FALSE, date='', years='five_year',range=5, update_run_date = FALSE,
                            withtag=FALSE, datatag='', use_share=FALSE, yl = year_list, region=NULL,
                            poly_ag = use_global_if_missing("poly_ag"),
                            zcol_ag = use_global_if_missing("zcol_ag")) {
  
  # Ensure str_match loaded
  str_match <- stringr::str_match
  
  if(withdate){
    if(date=='')
      rd=run_date
    if(date!='')
      rd=date
  } else {
    rd = run_date
  }
  
  # Load input data by indicator
  root <- ifelse(Sys.info()[1]=="Windows", <<<< FILEPATH REDACTED >>>>)
  if(use_share==FALSE) load_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
  if(use_share==TRUE) load_dir  <- <<<< FILEPATH REDACTED >>>>
  
  if(!withdate & !withtag) filename <- paste0(load_dir, indicator)
  if(withtag)              filename <- paste0(load_dir, indicator, datatag)
  if(withdate)             filename <- paste0(<<<< FILEPATH REDACTED >>>>)
  
  # try to see if an RDS exists, if so use that, if not use a csv
  if(file.exists(paste0(filename,'.RDS'))){
    message('READING INPUT DATA FROM RDS FILE')
    d <- readRDS(paste0(filename,'.RDS'))
  } else {
    message('READING INPUT DATA FROM CSV FILE')
    d <- fread(paste0(filename,'.csv'))
  }
  
  d$latitude  <- as.numeric(as.character(d$latitude))
  d$longitude <- as.numeric(as.character(d$longitude))
  message(nrow(d))
  
  # Remove odd long/lats for point data (when not using polygon resampling)
  if(!"point" %in% names(d)) {
    # assume only point data
    d$point = 1
  }
  
  d=d[d$latitude<=90 | (d$point == 0 & as.logical(poly_ag)),]
  d=d[d$latitude>=-90 | (d$point == 0 & as.logical(poly_ag)),]
  d=d[d$longitude<=180 | (d$point == 0 & as.logical(poly_ag)),]
  d=d[d$longitude>=-180 | (d$point == 0 & as.logical(poly_ag)),]
  d <- subset(d, !is.na(latitude) | (d$point == 0 & as.logical(poly_ag)))
  d <- subset(d, latitude!=0 | (d$point == 0 & as.logical(poly_ag)))
  message(nrow(d))
  
  # Check for necessary columns
  if(!(indicator %in% names(d))) stop(paste0("Your input data does not contain a column for your indicator: ", indicator))
  
  d <- as.data.table(d)
  
  # Change all "country" assignments to national level (in case subnational in the input data)
  if (nrow(d[grepl("[A-Z]*_[.]*", country),]) > 0) {
    subnat_countries <- unique(d[grepl("[A-Z]*_[.]*", country), country])
    warning(paste0("Changing subnational to national country codes for the following: ",
                   paste0(subnat_countries, collapse = ",")))
    d[grepl("[A-Z]*_[.]*", country), country := str_match(country,"([A-Z]*)_[.]*")[,2]]
  }
  
  
  # Subset to within modeling area
  if(!is.null(simple)){
    d$rowid <- 1:nrow(d)
    d$keep <- F
    dpoint <- d[d$point==1 | !as.logical(poly_ag), c("rowid", "longitude", "latitude")]
    coordinates(dpoint) <- c("longitude", "latitude")
    proj4string(dpoint) <- proj4string(simple)
    dpoint$keep <- !is.na(over(dpoint, as(simple, "SpatialPolygons")))
    d$keep[dpoint$rowid] <- dpoint$keep
  } else {
    d$keep <- T
  }
  
  # Keep polygon data that is in the region
  if(as.logical(poly_ag)) {
    if(!is.null(region)) {
      adm0_list <- get_adm0_codes(region, shapefile_version = modeling_shapefile_version)
      # get GAUL to iso mapping
      loc_codes <- get_location_code_mapping(shapefile_version=modeling_shapefile_version)
      loc_codes <- loc_codes[ADM_CODE %in% adm0_list, ]
      regs      <- loc_codes$ihme_lc_id
      d[country %in% regs & point==0, "keep"] <- T
    } else {
      warning("Missing region information. Will keep all polygon data")
      d[point==0, "keep"] <- T
    }
  }
  
  message(paste0(round(mean(d$keep), 2)*100, '% of input data in specified template'))
  d <- d[d$keep, ]
  
  if(agebin!=0)   d = d[age%in%agebin,]
  if(removeyemen) d = d[country!='Yemen' & country!='YEM',]
  
  # remap any years as needed
  if(years=='five_year') {
    d <- d[year >= 1998 & year <= 2002, year := 2000]
    d <- d[year >= 2003 & year <= 2007, year := 2005]
    d <- d[year >= 2008 & year <= 2012, year := 2010]
    d <- d[year >= 2013 & year <= 2017, year := 2015]
  }
  
  if (nrow(subset(d, year < min(yl))) > 0) {
    warning(paste0("Dropping all data before min(year_list) = ", min(yl), "..."))
    d <- subset(d, year >= min(yl))
  }
  if (nrow(subset(d, year > max(yl))) > 0) {
    warning(paste0("Dropping all data after max(year_list) = ", max(yl), "..."))
    d <- subset(d, year <= max(yl))
  }
  
  # add in a weight column if it does not exist
  if(!"weight" %in% names(d)) {
    warning("A 'weight' column does not exists so one is added with all 1s.")
    d[,weight := 1]
  }
  
  # creaste a weighted SS to base QTs on
  if(sum(c('N','weight') %in% colnames(d)) == 2) d[,weighted_n := N*weight]
  
  # check that there is some non aggregate data
  if(!is.null(zcol_ag)) {
    if(all( (d$point==0 & as.logical(poly_ag)) | !is.na(d[[zcol_ag]]))) {
      stop("There is only aggregate data. Currently the pipeline does not support modeling without some disaggregated data.")
    }
  } else {
    if(all(d$point == 0 & as.logical(poly_ag))) {
      stop("There is only aggregate data. Currently the pipeline does not support modeling without some disaggregated data.")
    }
  }
  
  # Save a copy
  if(update_run_date == TRUE) {
    if(dir.exists(paste0(<<<< FILEPATH REDACTED >>>>)) == TRUE) {
      existing_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
      new_try <- existing_dir
      index <- 0
      while(dir.exists(new_try)) {
        index <- index + 1
        new_try <- paste0(existing_dir, '_', index)
      }
      run_date <- paste0(run_date, '_', index)
      dir.create(new_try, showWarnings = FALSE)
      run_date_dir <- new_try
    }
    if(dir.exists(paste0(<<<< FILEPATH REDACTED >>>>)) == FALSE) {
      run_date_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
      dir.create(run_date_dir, showWarnings = FALSE)
    }
    write.csv(d, file=paste0(<<<< FILEPATH REDACTED >>>>))
    return(list(d, run_date))
  }
  
  if(update_run_date == FALSE) {
    if(agebin==0){
      run_date_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
      dir.create(run_date_dir, showWarnings = FALSE)
      write.csv(d, file=paste0(<<<< FILEPATH REDACTED >>>>))
    } else {
      run_date_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
      dir.create(run_date_dir, showWarnings = FALSE)
      write.csv(d, file=paste0(<<<< FILEPATH REDACTED >>>>))
    }
    return(d)
  }
  
}
