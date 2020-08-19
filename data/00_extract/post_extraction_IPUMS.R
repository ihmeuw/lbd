rm(list=ls())
cores <- 5
#Setup
folder_in <- '<<<< FILEPATH REDACTED >>>>'
folder_out <- '<<<< FILEPATH REDACTED >>>>'

today <- gsub("-", "_", Sys.Date())
#Load packages
packages <- c('haven', 'stringr', 'data.table', 'dplyr', 
              'magrittr', 'parallel', 'doParallel', 'feather')
packages <- lapply(packages, library, character.only=T)

message("Getting common column names")
#get the most recent pt and poly feathers and parse them for column names
pt <- '<<<< FILEPATH REDACTED >>>>' %>% 
  list.files(pattern="points", full.names=T) %>% 
  grep(pattern=".feather$", value=T)
pt <- pt[length(pt)] %>% feather_metadata
poly <- '<<<< FILEPATH REDACTED >>>>' %>% 
  list.files(pattern="points", full.names=T) %>% 
  grep(pattern=".feather$", value=T)
poly <- poly[length(poly)] %>% feather_metadata
pt_names <- pt[3] %>% unlist %>% names %>% gsub(pattern="types.", replacement="")
poly_names <- poly[3] %>% unlist %>% names %>% gsub(pattern="types.", replacement="")
noms <- c(pt_names, poly_names) %>% unique

message("List IPUMS dtas")
files <- list.files(folder_in, pattern="IPUMS_CENSUS", full.names=T)

message("Read in IPUMS Geo Codebook")
geo <- read.csv('<<<< FILEPATH REDACTED >>>>', stringsAsFactors = F, 
                encoding="windows-1252")
geo <- geo[, c("nid", "iso3", 'shapefile', 'location_code', 'lat', 
               'long', 'admin_level', 'point', 'geospatial_id', 
               'start_year', 'end_year')]

#define function to merge IPUMS files with geographies
ipums_merge <- function(dta, geo, folder_out, noms){
  dta <- read_dta(dta)
  #get dta survey info
  nid <- dta$nid[1] %>% as.character
  message(nid)
  iso3 <- dta$ihme_loc_id[1] %>% as.character
  year_start <- dta$year_start[1] %>% as.character
  year_end <- dta$year_end[1] %>% as.character
  survey_module <- dta$survey_module[1] %>% as.character
  
  #skip bad data
  has_pweight_not_hhweight <- "pweight" %in% names(dta) & !("hhweight" %in% names(dta))
  missing_all_weights <- !("pweight" %in% names(dta)) & !("hhweight" %in% names(dta))
  missing_wash <- !("t_type" %in% names(dta)) & 
    !("w_source_drink" %in% names(dta)) & !("sewage" %in% names(dta))
  if (has_pweight_not_hhweight){
    setnames(dta, "pweight", "hhweight")
  }
  if (missing_all_weights | missing_wash){
    return(nid)
  }
  
  m <- try(merge(dta, geo, by.x=c("nid", "ihme_loc_id", 'geospatial_id'), 
                 by.y=c("nid", 'iso3', 'geospatial_id'), all.x=T))
  
  if (class(m) == "try-error"){
    message(paste("Check try error", nid))
    return(nid)
  } else{
    outname <- paste("IPUMS_CENSUS", nid, survey_module, 
                     iso3, year_start, year_end, sep="_")
    m$survey_series <- m$survey_name
    m$iso3 <- m$ihme_loc_id
    orig_names <- names(m)
    new_names <- noms[!(noms %in% orig_names)]
    for (nam in new_names){
      m[, nam] <- NA
    }
    m <- m[, c("nid", 'survey_series', "year_start", "year_end", 
               "iso3", "lat", "long", "shapefile", "location_code", 
               "w_source_drink", "t_type", "sewage", "shared_san", 
               "mins_ws", "dist_ws", "hw_station", "hw_water", 
               "hw_soap", "hhweight", "hh_size", 'urban')]
    has_shp <- any(!is.na(m$shapefile)) & any(!is.na(m$location_code))
    has_lat_long <- any(!is.na(m$lat)) & any(!is.na(m$long))
    has_geography <- has_shp | has_lat_long
    if (!has_geography){
      message(paste("Check geography", nid, year_end))
      return(nid)
    }
    out <- paste0(folder_out, "/", outname, ".feather")
    write_feather(m, out)
    return(NULL)
  }
}

message("starting mclapply")
bad_nids <- mclapply(files, ipums_merge, geo=geo, folder_out=folder_out, 
                     noms=noms, mc.cores=cores) %>% unlist
write.csv(bad_nids, '<<<< FILEPATH REDACTED >>>>', row.names = F, na='')