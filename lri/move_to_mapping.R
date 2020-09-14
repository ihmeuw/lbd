#####################################################################
# Format LRI model outputs and move to mapping directory
#####################################################################
move_to_mapping <- function(map_date,
                            model_date,
                            measures,
                            include_aroc,
                            year_bounds,
                            stats,
                            admin_levels,
                            rasters,
                            include_counts){
  #(1) Setup ------------------------------------------------------------------
  #set paths
  path_to_inputs <- '<<<< FILEPATH REDACTED >>>>'
  
  
  for (measure in measures){
    
    map_dir <- '<<<< FILEPATH REDACTED >>>>'
    if (!dir.exists(map_dir)) dir.create(map_dir, recursive = TRUE)
    
    #(2) Rasters ---------------------------------------------------------------
    if (rasters){
      raster_dir <- paste0(path_to_inputs, '/pred_derivatives/global_rasters/')
      
      #move rasters
      for (stat in stats){
        move_command <- paste0('cp ', raster_dir, 'has_lri_', measure, '_', stat, '_raked_', year_bounds, '.tif ', map_dir, 'has_lri_', measure, '_', stat, '_raked_', year_bounds, '.tif')
        system(move_command)
      }
    }
    
    #(3) csvs ---------------------------------------------------------------
    for (admin in admin_levels){
      admin_code <- paste0('ADM', admin, '_CODE')
      agg_csv <- fread(paste0(path_to_inputs, '/pred_derivatives/admin_summaries/has_lri_admin_', admin, '_raked_', measure, '_summary.csv'))
      
      for (stat in stats){
        summary <- select(agg_csv, stat, year, admin_code) %>%
          rename(value = stat)
        write.csv(summary, paste0(map_dir, 'has_lri_', measure, '_', stat, '_raked_ad', admin, '.csv'))
      }
      
      if (include_aroc){
        aroc_csv <- fread(paste0(path_to_inputs, '/pred_derivatives/aroc/has_lri_admin_', admin, '_raked_', measure, '_aroc_summary.csv'))
        
        for (stat in stats){
          aroc <- select(aroc_csv, stat, year, admin_code) %>%
            rename(value = stat)
          write.csv(aroc, paste0(map_dir, 'has_lri_', measure, '_aroc_', stat, '_raked_ad', admin, '.csv'))
        }
        
        if (include_counts){
          count_csv <- fread(paste0(path_to_inputs, '/pred_derivatives/admin_summaries/has_lri_c_admin_', admin, '_raked_', measure, '_summary.csv'))
          
          for (stat in stats){
            count <- select(count_csv, stat, year, admin_code) %>%
              rename(value = stat)
            write.csv(count, paste0(map_dir, 'has_lri_', measure, '_', stat, '_counts_raked_ad', admin, '.csv'))
          }
        }
      }
    }
  }
}
