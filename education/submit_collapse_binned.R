package_list <- c("data.table", "gtools", "parallel", "SDMTools", "readstata13")


package_lib <- paste("<<<< FILEPATH REDACTED >>>>", paste(R.version$major, R.version$minor, sep="."), sep='/')
.libPaths(package_lib)
# Load your pacakges:
lapply(package_list, library,
       lib.loc        = package_lib,
       character.only = TRUE)


#Arguments
cores <- 20
update_only <- TRUE ## if TRUE, only crosswalk surveys that have been updated since update_since date
update_since <- as.POSIXct("2018-06-15") ## date after which any changed surveys should be crosswalked if update_only = T

##### Function to load edu data profiles and standardize names
## parameters: file = name of file located in prepped_data_profile_2 folder
load_edu_dta <- function(file) {
  file <- paste0("<<<< FILEPATH REDACTED >>>>",file)
  data <- data.table(read.dta13(file))
  if ("year_n" %in% names(data)) {
    data[, year := as.numeric(year_n)]
    data[, year_n := NULL]
  }
  data$file <- gsub("_profile", "", file)
  #data$file <- file
  underscores <- gregexpr("_", file)[[1]]
  data$nid <- as.numeric(substr(file, underscores[length(underscores)]+1,nchar(file)-4))
  print(file)
  return(data)
}

load_tab_dta <- function(file) {
  file <- paste0("<<<< FILEPATH REDACTED >>>>", file) 
  data <- data.table(read.csv(file))
  data$file <- file
  print(file)
  data <- unique(data[,c("bins", "type", "geospatial_id", "iso3", "level","ihme_loc_id","year","file","nid")])
  return(data)
}

## load all edu data profiles to figure out which surveys are binned and thus need to be crosswalked
data_dir <- "<<<< FILEPATH REDACTED >>>>"
files <- list.files(data_dir, pattern = ".dta")
data.list <- mclapply(files, load_edu_dta ,mc.cores=ifelse(Sys.info()[1]=="Windows", 1, cores))
edu.data <- rbindlist(data.list,fill=T)
edu.data <- edu.data[bins==1]


## load all tabulated data and construct profiles to figure out what needs to be crosswalked
tab.files <- list.files("<<<< FILEPATH REDACTED >>>>")
tab.files <- tab.files[tab.files != "metadata"]
tab.data <- rbindlist(mclapply(tab.files, load_tab_dta, mc.cores=ifelse(Sys.info()[1]=="Windows",1,cores)))
## drop any tab data we also have microdata for
tab.data <- tab.data[!tab.data$nid %in% edu.data$nid,]

edu.data <- rbind(edu.data,tab.data)

## only keep surveys post 1996 and in stages 1 or 2
edu.data <- edu.data[year > 1996]
stages <- read.csv("<<<< FILEPATH REDACTED >>>>")
st <- c("1", "2a", "2b")
stage_1_and_2_iso3 <- stages$iso3[stages$Stage %in% st]
edu.data <- edu.data[iso3 %in% stage_1_and_2_iso3]

## If we have update turned on, let's only re-do the files that have changed since update_since
if (update_only == TRUE) {
  files <- unique(edu.data$file)
  file_info <- file.info(files)
  file_info$file_name <- row.names(file_info)
  file_info <- as.data.table(file_info)
  updated_files <- file_info[mtime > update_since]
  edu.data <- edu.data[file %in% updated_files$file_name]
}

edu.data <- edu.data[level == 8| level ==7]

log_dir <- paste0("<<<< FILEPATH REDACTED >>>>",  format(Sys.Date(), "%Y_%m_%d"))
dir.create(log_dir, showWarnings = TRUE)

#launch one job per binned survey
for (i in 1:length(unique(edu.data$file))){
  dat <- edu.data[file == unique(edu.data$file)[i]]
  c.year <- unique(dat$year)
  c.level <- unique(dat$level)
  c.type <- unique(dat$type)
  cc.iso <- unique(dat$iso3)
  if (c.type=="BRFSS") {cc.iso <- "USA"}
  if (c.type=="UKGHS") {cc.iso <- "GBR"}
  c.file <- unique(dat$file)
  slashes <- gregexpr("/",c.file)[[1]]
  c.name <- substr(c.file, slashes[length(slashes)]+1, nchar(c.file)-4)

  system(paste0("qsub -o ", log_dir, " -e ", log_dir, " -N ", c.name  ," -P proj_geo_nodes -l geos_node=TRUE -pe multi_slot 10 <<<< FILEPATH REDACTED >>>>",
                 " ",c.year," ",c.level," " ,c.file))
  
}




