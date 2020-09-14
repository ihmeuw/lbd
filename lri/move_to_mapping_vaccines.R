#lri

map_date <- '2020_01_22'
model_date <- '2020_01_10_15_18_27'
measure <- 'mortality'
include_aroc <- FALSE
stats <- 'mean'
admin_levels <- 2
include_counts <- FALSE

#(1) Setup ------------------------------------------------------------------
#set paths
path_to_inputs <- '<<<< FILEPATH REDACTED >>>>'

  
  map_dir <- '<<<< FILEPATH REDACTED >>>>'
  if (!dir.exists(map_dir)) dir.create(map_dir, recursive = TRUE)
  
  
  #(3) csvs ---------------------------------------------------------------
  #Aggregation csvs: {indicator}_{measure}_{raked}_ad{x}.csv  with columns year, ADM{x}_CODE, and value
  for (admin in admin_levels){
    admin_code <- paste0('ADM', admin, '_CODE')
    agg_csv <- fread(paste0(path_to_inputs, '/pred_derivatives/admin_summaries/has_lri_lri_pneumo_admin_', admin, '_raked_', measure, '_summary.csv'))
    
    for (stat in stats){
      summary <- select(agg_csv, stat, year, admin_code) %>%
        rename(value = stat)
      write.csv(summary, paste0(map_dir, 'has_lri_pneumo_', measure, '_', stat, '_raked_ad', admin, '.csv'))
    }
  }

  #pcv cov
model_date <- '2019_12_17_16_03_24_raked'

#(1) Setup ------------------------------------------------------------------
#set paths
path_to_inputs <- '<<<< FILEPATH REDACTED >>>>'


map_dir <- '<<<< FILEPATH REDACTED >>>>'
if (!dir.exists(map_dir)) dir.create(map_dir, recursive = TRUE)


#(3) csvs ---------------------------------------------------------------
#Aggregation csvs: {indicator}_{measure}_{raked}_ad{x}.csv  with columns year, ADM{x}_CODE, and value
for (admin in admin_levels){
  admin_code <- paste0('ADM', admin, '_CODE')
  agg_csv <- fread('<<<< FILEPATH REDACTED >>>>')
  
  for (stat in stats){
    summary <- select(agg_csv, stat, year, admin_code) %>%
      rename(value = stat)
    write.csv(summary, paste0(map_dir, 'pcv3_cov_', stat, '_raked_ad', admin, '.csv'))
  }
  
  if (include_aroc){
    aroc_csv <- fread('<<<< FILEPATH REDACTED >>>>')
    
    for (stat in stats){
      aroc <- select(aroc_csv, stat, year, admin_code) %>%
        rename(value = stat)
      write.csv(aroc, paste0(map_dir, 'has_lri_', measure, '_aroc_', stat, '_raked_ad', admin, '.csv'))
    }
    
    if (include_counts){
      count_csv <- fread('<<<< FILEPATH REDACTED >>>>')
      
      for (stat in stats){
        count <- select(count_csv, stat, year, admin_code) %>%
          rename(value = stat)
        write.csv(count, paste0(map_dir, 'has_lri_', measure, '_', stat, '_counts_raked_ad', admin, '.csv'))
      }
    }
  }
}
}