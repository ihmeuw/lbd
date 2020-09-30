#########################################################################
# Read in malaria counterfactual csvs from Oxford, then merge and format
#########################################################################

# (1) Setup ---------------------------------------------------------------------
rm(list = ls())

ex_date <- '2019_09_18' #date of extraction in oxford_counterfactuals folder
run_date <- '2019_09_18' #date under which to save results in tridalys dir
scenarios <- c('act','all','irs','itn','raw')

#set directories
in_dir <- '<<< FILEPATH REDACTED >>>'
out_dir <- '<<< FILEPATH REDACTED >>>'
dir.create(out_dir, recursive = TRUE)

#read in link for shapefile version and subset to africa iso3s
link <- readRDS('<<< FILEPATH REDACTED >>>')
link <- link %>%
  select(ADM0_NAME, ADM0_CODE, ADM1_NAME, ADM1_CODE, ADM2_NAME, ADM2_CODE) %>%
  unique.data.frame()

# (2) Read in the files for one scenario ----------------------------------------
for (scen in scenarios){
  dir_to_search <- paste0(in_dir, scen)
  file_names <- list.files(dir_to_search, pattern = glob2rx('*infants*'))
  file_paths <- paste0(dir_to_search, '/', file_names)
  
  results <- list()
  for (n in 1:length(file_paths)){
    year_results <- fread(file_paths[n])
    results <- rbind(results, year_results)
  }
  
  results <- filter(results, Year %in% (2000:2017))
  results <- merge(results, link, by.x = c('Name','ISO3'), by.y = c('ADM2_NAME','ADM0_NAME')) %>%
    select(-'ID') %>%
    rename(c('ISO3' = 'ADM0_NAME', 'Name' = 'ADM2_NAME'))
  
  write.csv(results, file = paste0(out_dir, 'had_malaria_', scen, '_admin_2_raked_incidence_summary.csv'))
}

 
act <- fread('<<< FILEPATH REDACTED >>>')
all <- fread('<<< FILEPATH REDACTED >>>')
irs <- fread('<<< FILEPATH REDACTED >>>')
itn <- fread('<<< FILEPATH REDACTED >>>')
raw <- fread('<<< FILEPATH REDACTED >>>')
