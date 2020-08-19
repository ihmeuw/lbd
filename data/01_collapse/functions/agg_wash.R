agg_indi <- function(mydat = ptdat, var_family = indi_fam, 
                     dt_type = data_type, agg = agg_level,
                     condition = conditional) {
  
  if (var_family == 'water') {
      levels <- c('surface','imp','unimp',
                  'piped','piped_imp','piped_cw', 
                  'well_cw','well_imp','well_unimp',
                  'spring_cw','spring_imp','spring_unimp')
  }
  
  if (var_family == 'sani') {
    if (condition == 'conditional') {
      levels <- c('shared')
      mydat <- filter(mydat, imp == 1|latrine_imp == 1)
    } else {
      levels <- c('network','s_piped','imp','unimp','od',
                  'latrine_cw','latrine_imp','latrine_unimp',
                  'flush_imp','flush_unimp','flush_cw', 
                  'flush_imp_sewer','flush_imp_septic')
    }
    
  }
  mydat_0 <- mydat
  results <- list()
  for (i in levels) {
    mydat <- mydat_0

    message(paste("Aggregating",i))
    names(mydat)[which(names(mydat) == i)] <- 'indi'
    
    if (dt_type == 'pt') {
      if(agg == 'country') {
        
        mydatresults <- mydat %>% mutate(wt_indi = hhweight*indi*hh_size, 
                                         wt_denom = hhweight*hh_size,
                                         eff_n_num = hhweight*hh_size, 
                                         eff_n_denom = (hhweight^2)*hh_size, 
                                         year_binned = weighted) %>% 
          group_by(nid, iso3, survey_series, year_start) %>% 
          summarize(wtavg_indi = sum(wt_indi, na.rm = T)/sum(wt_denom, na.rm = T),
                    total_hh = ((sum(eff_n_num))^2)/sum(eff_n_denom)) %>%
          mutate(urban = NA)
        
      } else {
        
        mydatresults <- mydat %>% 
                          mutate(eff_n_num = hhweight*hh_size, 
                                 eff_n_denom = (hhweight^2)*hh_size) %>%
          
          group_by(id_short, nid, iso3, lat, long, survey_series, urban, 
                   year_start, int_year, year_median, shapefile, 
                   location_code) %>%
          
          
          summarize(wtavg_indi = (sum(indi*hh_size*hhweight, na.rm = TRUE)/
                                    sum(hh_size*hhweight, na.rm = T)),
                    
                    total_hh = ((sum(eff_n_num, na.rm = T))^2)/
                      sum(eff_n_denom, na.rm = T), 
                    
                    sum_of_sample_weights = sum(true_weight*hh_size, na.rm = T),
                    sum_old_N = sum(hh_size, na.rm = T))
      }
    }
    
    if (dt_type == 'poly') {
      
      if(agg == 'country') {
        
        mydatresults <- mydat %>% mutate(wt_indi = hhweight*indi*hh_size, 
                                         wt_denom = hhweight*hh_size,
                                         eff_n_num = hhweight*hh_size, 
                                         eff_n_denom = (hhweight^2)*hh_size) %>% 
          group_by(nid, iso3, survey_series, year_start) %>% 
          summarize(wtavg_indi = sum(wt_indi, na.rm = T)/sum(wt_denom, na.rm = T),
                    total_hh = ((sum(eff_n_num))^2)/sum(eff_n_denom)) %>%
          mutate(urban = NA)
        
      } else {
        mydatresults <- mydat %>% 
          mutate(# effective sample size by kish: OKAY
            eff_n_num = hhweight*hh_size, eff_n_denom = (hhweight^2)*hh_size) %>%
          
          group_by(id_short, nid, iso3, lat, long, survey_series, year_start, 
                   int_year, year_median, shapefile, location_code) %>%
          
          
          summarize(wtavg_indi = (sum(indi*hh_size*hhweight, na.rm = TRUE)/
                                    sum(hh_size*hhweight, na.rm = T)),
                    
                    total_hh = ((sum(eff_n_num, na.rm = T))^2)/
                      sum(eff_n_denom, na.rm = T), 
                    
                    sum_of_sample_weights = sum(true_weight*hh_size, na.rm = T),
                    sum_old_N = sum(hh_size, na.rm = T))%>%
          mutate(urban = NA)
      }
    } 
    
    names(mydatresults)[which(names(mydatresults) == 'wtavg_indi')] <- paste0(i)
    results[[length(results)+1]] <- mydatresults
    names(mydat)[which(names(mydat) == 'indi')] <- i
    
    
  }
  
  message("Merging all results...")
  mydat <- Reduce(function(x,y) merge(x,y,all = T),results)
  
  return(mydat)
}

