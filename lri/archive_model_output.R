#####################################################################
# Copy LRI run dates from given publication to archive folder on J
#####################################################################

archive_date <- '201909'
standard <- FALSE

if (standard){
  run_date_list <- c('2018_05_10_16_22_42',
                     '2018_05_10_16_22_42_covs',
                     '2018_05_10_16_22_42_gpcovs',
                     '2018_05_10_16_22_42_gponly',
                     '2018_05_10_16_22_42_stacking',
                     '2018_05_10_16_22_42_1000_draws',
                     '2018_05_10_15_58_20',
                     '2018_05_10_15_58_20_covs',
                     '2018_05_10_15_58_20_gpcovs',
                     '2018_05_10_15_58_20_gponly',
                     '2018_05_10_15_58_20_stacking',
                     '2018_05_10_15_58_20_1000_draws',
                     '2017_stage_1_updates')
  
  for (rd in run_date_list){
    print(paste0('Working on run date ', rd))
    cp_string <- '<<<< FILEPATH REDACTED >>>>'
    system(cp_string)
  }
} else {
  copy_from_list <- c('<<<< FILEPATH REDACTED >>>>','<<<< FILEPATH REDACTED >>>>')
  copy_to_list <- c('<<<< FILEPATH REDACTED >>>>',
                    '<<<< FILEPATH REDACTED >>>>')

  for (n in 1:length(copy_from_list)){
    print(paste0('Working on directory ', copy_from_list[n]))
    cp_string <- paste0('cp -r ', copy_from_list[n], ' ', copy_to_list[n])
    system(cp_string)
  }
}
