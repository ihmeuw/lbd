#####################################################################
# Copy LRI run dates from given publication to archive folder on J
#####################################################################
archive_date <- '201909'

run_dates_to_find <- c('2018_05_10_16_22_42',
                       '2018_05_10_15_58_20',
                       '2017_stage_1_updates')

dir_to_search <- '<<<< FILEPATH REDACTED >>>>'
save_dir <- '<<<< FILEPATH REDACTED >>>>'
dir.create(save_dir)

run_date <- run_dates_to_find[1]
file_list <- list.files(dir_to_search, pattern = run_date)

for (run_date in run_dates_to_find){
  print(paste0('Working on run date ', run_date))
  file_list <- list.files(dir_to_search, pattern = run_date)
  for (file in file_list){
    cp_string <- paste0('cp -r ', dir_to_search, '/', file, ' ', save_dir, file)
    system(cp_string)
  }
}