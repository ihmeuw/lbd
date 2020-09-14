# Compare model results between two run dates by admin

adm_compare <- function(indicator,
                        indicator_group,
                        rake_list,
                        admin_lvls,
                        measure_list,
                        age,
                        holdout,
                        run_date_1,
                        run_date_2,
                        region_list,
                        collapse_years = F) {
  #define outdir
  out_dir <- '<<<< FILEPATH REDACTED >>>>'
  
  #load packages
  library(ggplot2)
  library(dplyr)
  
  #loop over various comparisons and create graphs
  for (measure in measure_list){
    for (admin_lvl in admin_lvls){
      for (raked in rake_list){
        #load in RData of admin-level means
        admin_old <- read.csv('<<<< FILEPATH REDACTED >>>>')
        
        admin_new <- read.csv('<<<< FILEPATH REDACTED >>>>')
        
        #change column names to prepare for merge
        admin_new$mean_new <- admin_new$mean
        admin_old$mean_old <- admin_old$mean
        admin_old$mean <- NULL
        admin_new$mean <- NULL
        
        admin_new$lower_new <- admin_new$lower
        admin_old$lower_old <- admin_old$lower
        admin_old$lower <- NULL
        admin_new$lower <- NULL
        
        admin_new$upper_new <- admin_new$upper
        admin_old$upper_old <- admin_old$upper
        admin_old$upper <- NULL
        admin_new$upper <- NULL
        
        admin_new$cirange_new <- admin_new$cirange
        admin_old$cirange_old <- admin_old$cirange
        admin_old$cirange <- NULL
        admin_new$cirange <- NULL
        
        #merge old and new admin summaries by ADM0_CODE, ADM0_NAME, year
        if (admin_lvl == 0){
          admin_compare <- merge(admin_old, admin_new, by = c("ADM0_CODE", "ADM0_NAME", "region", "year"), all.x = T, all.y = T)
        } else if (admin_lvl == 1){
          admin_compare <- merge(admin_old, admin_new, by = c("ADM0_CODE", "ADM0_NAME", "ADM1_CODE", "ADM1_NAME", "region", "year"), all.x = T, all.y = T)
        } else if (admin_lvl == 2){
          admin_compare <- merge(admin_old, admin_new, by = c("ADM0_CODE", "ADM0_NAME", "ADM1_CODE", "ADM1_NAME", "ADM2_CODE", "ADM2_NAME", "region", "year"), all.x = T, all.y = T)
        }
        
        if (!dir.exists(out_dir)) {
          dir.create(out_dir)
        }
        
      if (collapse_years == F){
        pdf(paste0(out_dir, '/compare_', ifelse(raked, "raked_", "unraked_"), measure, '_adm_', admin_lvl, '_', run_date_1, '_to_', run_date_2, '.pdf'), width = 12, height = 10)
      } else if (collapse_years == T){
        pdf(paste0(out_dir, '/compare_all_years_', ifelse(raked, "raked_", "unraked_"), measure, '_adm_', admin_lvl, '_', run_date_1, '_to_', run_date_2, '.pdf'), width = 12, height = 10)
      }
        
        for (reg in region_list){
          
          #ggplot a scatterplot of mean_new vs. mean_old, facetted by year
          if (collapse_years == F){
            admin_compare_years <- filter(admin_compare, year %in% seq(2000, 2015, by = 5) & region == reg)
            graph <- ggplot(data = admin_compare_years, aes(mean_old, mean_new)) + geom_point(aes(color = ADM0_NAME)) + 
              facet_wrap(~year) + geom_abline(intercept = 0, slope = 1) + ggtitle(paste0(indicator, ' ', measure, ifelse(raked, ' raked', ' unraked'), ': ', reg, ' ADM ', admin_lvl))
          } else if (collapse_years == T){
            admin_compare_years <- filter(admin_compare, region == reg)
            graph <- ggplot(data = admin_compare_years, aes(mean_old, mean_new)) + geom_point(aes(color = ADM0_NAME)) + 
               geom_abline(intercept = 0, slope = 1) + ggtitle(paste0(indicator, ' ', measure, ifelse(raked, ' raked', ' unraked'), ': ', reg, ' ADM ', admin_lvl))
          }
        
          plot(graph)
        }
        dev.off()
      }
    }
  }
}

admin_correlation_coefs <-function(indicator,
                                   indicator_group,
                                   raked,
                                   admin_lvl,
                                   measure,
                                   age,
                                   holdout,
                                   run_date_1,
                                   run_date_2,
                                   region_list,
                                   collapse_years = F){
  #load in RData of admin-level means
  admin_old <- read.csv('<<<< FILEPATH REDACTED >>>>')
  
  admin_new <- read.csv('<<<< FILEPATH REDACTED >>>>')
  
  #change column names to prepare for merge
  admin_new$mean_new <- admin_new$mean
  admin_old$mean_old <- admin_old$mean
  admin_old$mean <- NULL
  admin_new$mean <- NULL
  
  admin_new$lower_new <- admin_new$lower
  admin_old$lower_old <- admin_old$lower
  admin_old$lower <- NULL
  admin_new$lower <- NULL
  
  admin_new$upper_new <- admin_new$upper
  admin_old$upper_old <- admin_old$upper
  admin_old$upper <- NULL
  admin_new$upper <- NULL
  
  admin_new$cirange_new <- admin_new$cirange
  admin_old$cirange_old <- admin_old$cirange
  admin_old$cirange <- NULL
  admin_new$cirange <- NULL
  
  #merge old and new admin summaries by ADM0_CODE, ADM0_NAME, year
  if (admin_lvl == 0){
    admin_compare <- merge(admin_old, admin_new, by = c("ADM0_CODE", "ADM0_NAME", "region", "year"), all.x = T, all.y = T)
  } else if (admin_lvl == 1){
    admin_compare <- merge(admin_old, admin_new, by = c("ADM0_CODE", "ADM0_NAME", "ADM1_CODE", "ADM1_NAME", "region", "year"), all.x = T, all.y = T)
  } else if (admin_lvl == 2){
    admin_compare <- merge(admin_old, admin_new, by = c("ADM0_CODE", "ADM0_NAME", "ADM1_CODE", "ADM1_NAME", "ADM2_CODE", "ADM2_NAME", "region", "year"), all.x = T, all.y = T)
  }
  
  results_list <- data_frame()
  
  for (n in 1:length(region_list)){
    reg <- region_list[n]
 
    if (collapse_years == F){
      for (yr in year_list){
        admin_compare_years <- filter(admin_compare, region == reg & year == yr)
        cor <- cor(admin_compare_years$mean_new, admin_compare_years$mean_old, use = 'complete.obs')
        results <- data_frame('reg' = reg, 'year' = yr, 'cor' = cor)
        results_list <- rbind(results_list, results)
      }
      
    } else if (collapse_years == T){
      admin_compare_years <- filter(admin_compare, region == reg)
      cor <- cor(admin_compare_years$mean_new, admin_compare_years$mean_old, use = 'complete.obs')
      results <- data_frame('reg' = reg, 'cor' = cor)
      results_list <- rbind(results_list, results)
    }
  }
  return(results_list)
}

