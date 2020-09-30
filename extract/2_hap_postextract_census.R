rm(list=ls())

#Define values
#Redacted

#Setup
#Redacted

#Load packages
package_lib    <- file.path(h, '_code/_lib/pkg')
  .libPaths(package_lib)
pacman::p_load(data.table, dplyr, haven, feather, fst, googledrive, magrittr, parallel, doParallel, readxl, stringr)

#read stage document
stages <- file.path('#Redacted/stage_master_list.csv') %>% fread #read info about stages

#timestamp
today <- Sys.Date() %>% gsub("-", "_", .)

message("Getting common column names")
if (topic == "hap" & geos){
  #get the most recent pt and poly files and parse them for column names
  pt_names <- paste0(l, "#Redacted") %>% 
    list.files(pattern="points", full.names=T) %>% 
    grep(pattern=".fst$", value=T) %>% 
    .[length(.)] %>% 
    read_fst(., from = 1, to = 2) %>% 
    names
  poly_names <- paste0(l, "#Redacted") %>% 
    list.files(pattern="poly", full.names=T) %>% 
    grep(pattern=".fst$", value=T) %>% 
    .[length(.)] %>% 
    read_fst(., from = 1, to = 2) %>% 
    names
  noms <- c(pt_names, poly_names) %>% unique
} else{
  #Redacted
}

message("List IPUMS dtas")
extractions <- list.files(folder_in, pattern="IPUMS", full.names=T)

#Change to handle batch extractions by only reading in those IDs that have been extracted by #Redacted
#Redacted
codebook <- read_xlsx('hap.xlsx', sheet='codebook') %>% as.data.table
#create output name, note that we need to remove the leading info on some of the survey names(take only str after /)
codebook[, output_name := paste0(ihme_loc_id, '_', tools::file_path_sans_ext(basename(survey_name)), '_', year_start, '_', year_end, '_', nid, '.csv')]
#subset the codebook to ONLY the files that our extractors have worked on
codebook <- codebook[assigned %in% extractor_ids]
#get the names of all the new files
new.files <- codebook[, output_name] %>% unique %>% paste(., collapse="|")
extractions <- grep(new.files, extractions, invert=F, value=T)

#define function to merge IPUMS files with geographies
ipums_merge <- function(file, geo, folder_out, noms){
  
  dt <- fread(file)

  #get survey info
  nid <- dt$nid[1] %>% as.character
  message(nid)
  iso3 <- dt$ihme_loc_id[1] %>% as.character
  year_start <- dt$year_start[1] %>% as.character
  year_end <- dt$year_end[1] %>% as.character
  survey_module <- dt$survey_module[1] %>% as.character
  
  #use new geocodebook database 
  source("#Redacted/geocodebook_functions.R")
  geo <- get_geocodebooks(nids = nid)
  
  #skip bad data
  has_pweight_not_hhweight <- "pweight" %in% names(dt) & !("hhweight" %in% names(dt))
  missing_all_weights <- !("pweight" %in% names(dt)) & !("hhweight" %in% names(dt))
  missing_gid <- !("geospatial_id" %in% names(dt))
  missing_hap <- !("cooking_fuel" %in% names(dt)) & !(grepl('housing', names(dt)) %>% any)
  
  if (has_pweight_not_hhweight) setnames(dt, "pweight", "hhweight")

  if (missing_all_weights | missing_hap | missing_gid) return(nid)
  
  if (!missing_gid) dt[, geospatial_id := as.character(geospatial_id)] #force geospatial IDs to character to match the sheet
  m <- try(merge(dt, geo, by.x=c("nid", "ihme_loc_id", 'geospatial_id'), by.y=c("nid", 'iso3', 'geospatial_id'), all.x=T))

  if (class(m) == "try-error") {message(paste("Check try error", nid)); return(nid)}
  else{
    outname <- paste("IPUMS_CENSUS", nid, survey_module, iso3, year_start, year_end, sep="_")
    m$survey_series <- m$survey_name
    m$iso3 <- m$ihme_loc_id
    orig_names <- names(m)
    new_names <- noms[!(noms %in% orig_names)]
    for (nam in new_names){
      m[, nam] <- NA
    }
    has_shp <- any(!is.na(m$shapefile)) & any(!is.na(m$location_code))
    has_lat_long <- any(!is.na(m$lat)) & any(!is.na(m$long))
    has_geography <- has_shp | has_lat_long
    if (!has_geography){message(paste("Check geography", nid, year_end)); return(nid)}

    out <- paste0(folder_out, "/", outname, ".fst")
    write.fst(m, out)
    return(NULL)
  }
}

message("starting mclapply")
bad_nids <- mclapply(extractions, ipums_merge, geo=geo, folder_out=folder_out, noms=noms, mc.cores=cores) %>% unlist
write.csv(bad_nids, paste0(folder_out, "/fix_these_nids.csv"), row.names = F, na='')


#####################################################################
#######Find broken extractions###################
#####################################################################

files <- list.files('#Redacted', pattern = 'IPUMS_CENSUS')
files <- substr(files, 1, nchar(files)-4)
files <- strsplit(files,'_')
files[6] #TODO jank
for(f in 1:length(files)){
  temp <- unlist(files[f])
  temp <- temp[6]
  files[f] <- temp
}
files <- unlist(files)
#find nids in codebook not in this list ^
codebook.nids <- codebook[!(year_end < 2000 | ihme_loc_id %in% stages[Stage==3, iso3]), nid] %>% unique
bad_ipums <- codebook.nids[!(codebook.nids %in% files)]

#check against bad nids csv
broken_ipums <- bad_nids[!(bad_nids %in% bad_ipums)]

##Ipums failed extractions
write.csv(broken_ipums, paste0(folder_out, '/failed_extractions_ipums.csv'))