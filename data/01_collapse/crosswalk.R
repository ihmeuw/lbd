rm(list = ls())

repo <- '<<<< FILEPATH REDACTED >>>>'

#get latest file
files <- file.info(list.files('<<<< FILEPATH REDACTED >>>>', 
                              pattern = '*.feather', full.names=TRUE))
files <- files[with(files, order(as.POSIXct(ctime), decreasing = TRUE)), ]
latest_postextraction <- unlist(strsplit(rownames(files)[1], "/"))
latest_postextraction <- latest_postextraction[length(latest_postextraction)]
input_version <- gsub('.feather', '', latest_postextraction)
input_version <- substr(input_version,(nchar(input_version) - 9),
                        nchar(input_version))


setwd(repo)
source('functions/cw_indi.R') #cw functions
#generats sb2
source('functions/tabulated_data_fix.R') #handle India tabulated data

library(dplyr)
library(feather)
library(data.table)

#grab data
setwd('<<<< FILEPATH REDACTED >>>>')
points <- read_feather(paste0('ptdat_sani_unconditional__', 
                              input_version,'.feather'))
poly <- read_feather(paste0('polydat_sani_unconditional__', 
                            input_version,'.feather'))
#Grab data in table form not put through pipeline
tabs <- fread('tabulated_data/sani.csv')

setwd(paste0('<<<< FILEPATH REDACTED >>>>', input_version))
ipums <- list.files(pattern = 'sani_')
ipums <- lapply(ipums, read_feather)
ipums <- do.call(rbind, ipums)

ipums$location_code <- as.character(ipums$location_code)
alldat <- as.data.frame(bind_rows(points, poly, ipums))
alldat <- as.data.frame(bind_rows(alldat, tabs))
alldat$iso3 <- substr(alldat$iso3, 1, 3)
cw_dat <- cw_sani(alldat)
today <- gsub("-", "_", Sys.Date())

cw_dat <- subset(cw_dat, year_start > 1999)
sb2 <- rename(sb2, imp = s_imp_other, piped = s_piped, network = s_network, 
              unimp = s_unimp, od = s_od,
              year_start = year, sum_old_N = total_hh) %>%
  dplyr::select(-imp_cw)

write_feather(cw_dat, 
			  paste0('<<<< FILEPATH REDACTED >>>>',
			  	     today, '.feather'))

###
setwd('<<<< FILEPATH REDACTED >>>>')
points <- read_feather(paste0('ptdat_water_unconditional__',
                              input_version,'.feather'))
poly <- read_feather(paste0('polydat_water_unconditional__',
                            input_version,'.feather'))
alldat <- as.data.frame(bind_rows(points, poly))
tabs <- fread('tabulated_data/water.csv')

alldat <- as.data.frame(bind_rows(alldat, tabs))
alldat$iso3 <- substr(alldat$iso3, 1, 3)
cw_dat <- cw_water(alldat)

cw_dat <- subset(cw_dat, year_start > 1999)

write_feather(cw_dat, 
			  paste0('<<<< FILEPATH REDACTED >>>>',
			  	     today, '.feather'))

