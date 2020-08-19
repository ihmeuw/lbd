rm(list = ls())
# Set library and load packages
## drive locations
commondir      <- sprintf('<<<< FILEPATH REDACTED >>>>')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))


# TBD: Remve all 'setwd()'
core_repo <- repo <-  '<<<< FILEPATH REDACTED >>>>'
setwd(repo)

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = repo)
library(feather)

proj <- '-P proj_geospatial'
user <- Sys.info()[['user']]

setwd('<<<< FILEPATH REDACTED >>>>')
indicators <- c('sani','water')
run_date <- Sys.Date()

for (indic in indicators) {
  if (indic == 'water') {
    polydat <- read_feather('<<<< FILEPATH REDACTED >>>>')
  } else {
    polydat <- read_feather('<<<< FILEPATH REDACTED >>>>')
  }
  
  polydat <- subset(polydat, is.na(lat) & !is.na(shapefile) & !is.na(location_code))
  for (shp in unique(polydat$shapefile)) { 
    jname <- paste(indic, shp, sep = "_")
    mythreads <- '1'
    mymem <- '5G'
    sys.sub <- paste0("qsub ",proj,paste0(" -e <<<< FILEPATH REDACTED >>>> -o <<<< FILEPATH REDACTED >>>> "),
                      "-cwd -N ", jname, " ", "-l fthread=", mythreads, " ", "-l m_mem_free=", mymem, 
                      ' -q long.q -l h_rt=06:00:00 -l archive=TRUE')
    script <- "child.R"
    r_shell <- '<<<< FILEPATH REDACTED >>>>'
    
    args <- paste(shp, indic, run_date)
    system(paste(sys.sub, r_shell, script, args)) 
  }
}
