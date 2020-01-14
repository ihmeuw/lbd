#### Data Validation Script ####
#### Aniruddha Deshpande, adesh@uw.edu ####

# Clear environment of all objects
rm(list = ls())

#### User Defined Inputs ####
# Define path to save logs to
log_path <- 'C:/Users/adesh/Desktop/parsing.txt'

# Define filepath to save country data plots to
plot_path <- 'C:/Users/adesh/Desktop/country_data_plots.pdf'

# Indicator name
indi_name <- 'w_piped'

# Drive where input data is stored ('j' or 'share')
data_drive <- 'j'

# Variable names. Dataset should have the following variables:
# 1. 'country': indicates a country via ISO3 code
# 2. 'latitude': latitude of associated data record
# 3. 'longitude': longitude of associated data record
# 4. 'year': indicates the year associated with the data record
# 5. 'N': indicates sample size
# 6. indicator: Should be the same as the 'indi_name' supplied above
# 7. 'prop': Should be the indicator expressed as a proportion
# if any variables aren't named as above please supply the names below

country_var_name <- 'iso3'
latitude_var_name <- 'latitude'
longitude_var_name <- 'longitude'
year_var_name <- 'year'
N_var_name <- 'N'
indicator_var_name <- indi_name

# Do you have a variable that expresses the indicator as a proportion?
prop_present <- 'T'
prop_var_name <- 'prop'

#### Set Up Environment ####
# Detach all packages
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices",
                      "package:utils","package:datasets","package:methods",
                      "package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) 
    detach(package, character.only=TRUE)
}
detachAllPackages()

# Set library path
if(Sys.info()[1]!="Windows") {
  root <- "/home/j/"
  package_lib <- ifelse(grepl("geos", Sys.info()[4]),
                        paste0(root,'temp/geospatial/geos_packages'),
                        paste0(root,'temp/geospatial/packages'))
  .libPaths(package_lib)
} else {
  package_lib <- .libPaths()
  root <- 'J:/'
}

message(paste0('Loading packages from ',package_lib))

# Load Packages
packages <- c('dplyr', 'readr', 'tidyr', 'ggplot2', 'tmap', 'sp', 'RColorBrewer', 'gridExtra')
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(packages, library, character.only = T)

# Load Data & Write Data Parsing Log
parse_log <- file(log_path, open = 'w')
sink(parse_log)
sink(parse_log, type = 'message', append = T)
message('1. Data Parsing')
mydat <- read_csv(ifelse(data_drive == 'j',
                         paste0(root,'WORK/11_geospatial/10_mbg/input_data/',indi_name,'.csv'),
                         paste0('/share/geospatial/mbg/input_data/',indi_name,'.csv')
        ))
names(mydat)[which(names(mydat) == country_var_name)] <- 'country'
names(mydat)[which(names(mydat) == latitude_var_name)] <- 'latitude'
names(mydat)[which(names(mydat) == longitude_var_name)] <- 'longitude'
names(mydat)[which(names(mydat) == year_var_name)] <- 'year'
names(mydat)[which(names(mydat) == N_var_name)] <- 'N'
names(mydat)[which(names(mydat) == indi_name)] <- 'indi'

if (prop_present) {
  names(mydat[which(names(mydat) == prop_var_name)]) <- 'prop'
} else {
  mydat <- mutate(mydat, prop = indi/N)
}

#### Data Validation ####
### 1. Space ###
message('\n 2.1 Space Records Validation')

# Checking validity of space records
inv_lat <- length(mydat$latitude[which(mydat$latitude < -90 &
                                        mydat$latitude > 90)])
inv_long <- length(mydat$longitude[which(mydat$longitude < -180 &
                                         mydat$longitude > 180)])

# Defining standard list of iso3 codes from UN website
iso_list <- c('AFG','ALB','DZA','AND','AGO','ATG','ARG','ARM','AUS','AUT','AZE','BHS',
               'BHR','BGD','BRB','BLR','BEL','BLZ','BEN','BTN','BOL','BIH','BWA','BRA',
               'BRN','BGR','BFA','BDI','CPV','KHM','CMR','CAN','CAF','TCD','CHL','CHN',
               'COL','COM','COG','COK','CRI','HRV','CUB','CYP','CZE','CIV','PRK','COD',
               'DNK','DJI','DMA','DOM','ECU','EGY','SLV','GNQ','ERI','EST','ETH','FRO',
               'FJI','FIN','FRA','GAB','GMB','GEO','DEU','GHA','GRC','GRD','GTM','GIN',
               'GNB','GUY','HTI','HND','HUN','ISL','IND','IDN','IRN','IRQ','IRL','ISR',
               'ITA','JAM','JPN','JOR','KAZ','KEN','KIR','KWT','KGZ','LAO','LVA','LBN',
               'LSO','LBR','LBY','LTU','LUX','MDG','MWI','MYS','MDV','MLI','MLT','MHL',
               'MRT','MUS','MEX','FSM','MCO','MNG','MNE','MAR','MOZ','MMR','NAM','NRU',
               'NPL','NLD','NZL','NIC','NER','NGA','NIU','NOR','OMN','PAK','PLW','PAN',
               'PNG','PRY','PER','PHL','POL','PRT','QAT','KOR','MDA','ROU','RUS','RWA',
               'KNA','LCA','VCT','WSM','SMR','STP','SAU','SEN','SRB','SYC','SLE','SGP',
               'SVK','SVN','SLB','SOM','ZAF','SSD','ESP','LKA','SDN','SUR','SWZ','SWE',
               'CHE','SYR','TJK','THA','MKD','TLS','TGO','TKL','TON','TTO','TUN','TUR',
               'TKM','TUV','UGA','UKR','ARE','GBR','TZA','USA','URY','UZB','VUT','VEN',
               'VNM','YEM','ZMB','ZWE')

# Identifying non-standard iso3 codes and matching to standard ones
prob_iso <- unique(mydat$country[which(!(mydat$country %in% iso_list))])
for (i in prob_iso) {
  rep_iso <- iso_list[grep(substr(prob_iso, 1, 3), iso_list)]
  mydat$country[which(mydat$country == i)] <- rep_iso
}

message(paste(
  'There are', inv_lat, 'records with invalid latitudes.',
  '\n There are', inv_long, 'records with invalid longitudes.',
  '\n Country plots of data are being stored at: \n', plot_path,
  '\n Verify all points fall within country borders.',
  '\n There are', length(prob_iso), 'non-standard ISO3 codes.'
  ))
if(length(prob_iso)>0) {
  message(paste('They are', prob_iso, 
                'They are replaced using closest matches in list.'))
  }

### 2. Time ###
message('\n 2.2 Time Records Validation')

# Checking if all time records are integers
year_int <- F %in% (mydat$year == as.integer(mydat$year))

if (year_int) {
  message('Not all time records are integers.')
} else {message('All time records are integers.')}

# Checking validity of time records
old_thresh <- 1990
old_years <- length(mydat$year[which(mydat$year < old_thresh)])

neg_years <- length(mydat$year[which(mydat$year < 0)])

message(paste(
  'There are', old_years, 'records before', old_thresh, '.',
  'There are', neg_years, 'records with negative time.'
))

# Checking association of nid with year
nid_year <- count(mydat, country, year, nid)
nid_reps <- as.character(unique(nid_year[duplicated(nid_year$nid),'nid']))
time_probs <- filter(nid_year, nid %in% nid_reps)
if (length(unique(time_probs$nid)) > 0) {
  message(paste(
      'There are', length(unique(time_probs$nid)),
      'NIDs associated with multiple years. They are located in:',
      unique(time_probs$country),'.'
  ))
  } else {message('All NIDs are associated with a single year.')}

### 3. Denominator ###
message('\n 2.3 Denominator Records Validation')

# Checking if all records are integers
n_int <- F %in% (mydat$N == as.integer(mydat$N))

if (n_int) {
  message('Not all denominator records are integers.')
} else {message('All denominator records are integers.')}

# Checking validity of denominator records
neg_denom <- length(mydat$N[which(mydat$N < 1)])
message(paste(
  'There are', neg_denom, 
  'denominator records with 0 or negative values.'
))

### 4. Numerator ###
message('\n 2.4 Successes Records Validation')

# Checking if all records are integers
ind_int <- F %in% (mydat$N == as.integer(mydat$indi))
if (ind_int) {
  message('Not all successes records are integers.')
} else {message('All successes records are integers.')}

# Checking validity of success records
neg_ind <- length(mydat$indi[which(mydat$indi < 0)])
message(paste(
  'There are', neg_ind, 
  'success records with negative values.'
))

inv_ind <- length(which(mydat$indi > mydat$N))
message(paste(
  'There are', inv_ind,
  'success records with values greater than the respective denominator.'
))

sink(NULL)
sink(NULL, type = 'message')


### Create country data plots ###
# Write plots to PDF
pdf(plot_path)
for (i in unique(mydat$country)) {
  # Subset dataset to country level data
  dxdat <- mydat[which(mydat$country == i),]

  # Write summary statistics for numerator and denominator
  plot.new()
  print(grid.table(data.frame(stat = names(summary(dxdat$N)), 
                                           N = as.numeric(summary(dxdat$N)))))
  plot.new()
  print(grid.table(data.frame(stat = names(summary(dxdat$prop)), 
                                           prop = as.numeric(summary(dxdat$prop)))))

  # Write total number of people represented by each survey
  print(grid.table(dxdat %>% mutate(new_N = weight*N) %>% 
    group_by(country, year, survey_series, nid) %>% summarize(total_people = sum(new_N))))

  # Make histograms of numerator and denominator
  print(hist(dxdat$N))
  print(hist(dxdat$prop))

  # Make a scatter plot of numerator vs denominator with coloring by survey_series
  print(ggplot(dxdat) + geom_point(aes(x = N, y = prop, color = survey_series)) + 
    ggtitle(i) + theme_bw())

  # Plot the data against the shapefile and color coded for survey_series
  data(World)
  dxdat_sp <- SpatialPointsDataFrame(coords = data.frame(dxdat$longitude, dxdat$latitude),
                                   data = dxdat)
  print(tm_shape(World[which(World$iso_a3 == i),]) + tm_polygons() + tm_legend(title = i) +
    tm_shape(dxdat_sp) + tm_dots(col = 'survey_series', palette = 'Set1'))
}
dev.off()
