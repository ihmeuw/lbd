# Function to see data changes
library(tidyr)
library(dplyr)

geomatched_data_changes <- function(old_data, new_data) {
  cat("\n\n*****NIDs that have changed:\n")
  for (nn in sort(intersect(old_data$nid, new_data$nid))) {
    temp_old <- old_data[nid == nn, mget(c("nid", "country", "source", "year", "geospatial_id",
                                      "point","latitude","longitude","shapefile","indicator","pweight", 
                                      "age_year", "hh_id", "line_id", "location_code", "psu"))] %>%
      arrange(nid, country, year, geospatial_id, point, latitude, longitude, indicator, pweight, age_year, hh_id, line_id, location_code, psu)
    
    temp_new <- new_data[nid == nn, mget(c("nid", "country", "source", "year", "geospatial_id",
                                           "point","latitude","longitude","shapefile","indicator","pweight", 
                                           "age_year", "hh_id", "line_id", "location_code", "psu"))] %>%
      arrange(nid, country , year, geospatial_id, point, latitude, longitude, indicator, pweight, age_year, hh_id, line_id, location_code, psu)
    
    # If not equal, run some data checks
    if (!are_equal(temp_old, temp_new)) {
      cat(paste("\n\n", nn, "- data changed\n"))
      
      # Sample size and prevalence of new and old nids, make sure it is identical
      old_nid <-
        temp_old %>%
        dplyr::summarize(ind_prev = round(100 * weighted.mean(indicator, pweight), 3),
                  sample_size = n())
    
      new_nid <-
        temp_new %>%
        dplyr::summarize(ind_prev = round(100 * weighted.mean(indicator, pweight), 3),
                  sample_size = n())
      
      if (old_nid$ind_prev != new_nid$ind_prev) {
        cat(paste0("\n\n !! Weighted prevalence changed from ", old_nid$ind_prev, "% to ", new_nid$ind_prev, "%, this should not happen"))
      } else {
        cat("\n\nPrevalence remained the same")
      }
      
      if (old_nid$sample_size != new_nid$sample_size) {
        cat(paste0("\n\n !! Sample size changed from ", old_nid$sample_size, " (old data) to ", new_nid$sample_size, "; (new data) a difference of ", 
                   new_nid$sample_size - old_nid$sample_size, " people (", round(100 * abs(new_nid$sample_size - old_nid$sample_size) / nrow(temp_new), 1), "%)",
                   "\nthis should not happen"))
      } else {
        cat("\n\nSample size remained the same")
      }
      
      # Points
      old_points <- data.table(point = c(0, 1), n = c(nrow(filter(temp_old, point == 0)), nrow(filter(temp_old, point == 1))))
      new_points <- data.table(point = c(0, 1), n = c(nrow(filter(temp_new, point == 0)), nrow(filter(temp_new, point == 1))))
      
      if (are_equal(old_points, new_points)) {
        cat(paste0("\n\nGeographic information remained the same: ", new_points[point == 1, n], " people mapped to points and ", new_points[point == 0, n], " mapped to polygons\n\n"))
      } else {
        cat(paste0("\n\n !! Geographic information changed:\nthere were ", old_points[point == 1, n], " people mapped to points and ", old_points[point == 0, n],
                   " mapped to polygons in the old data\nand now there are ", new_points[point == 1, n], " people mapped to points and ", new_points[point == 0, n], " mapped to polyons in the new data \n\n"))
      }
      
      # Additional checks if lengths are the same but still differences
      if (are_equal(dim(temp_old), dim(temp_new))){
        cat("all.equal function run for additional diagonstic information (current = new data, target = old data)\n\n")
        print(all.equal(temp_old, temp_new))
      }
      
      # are_equal(temp_old$location_code, temp_new$location_code)
      # are_equal(temp_old$age_year, temp_new$age_year)
      # setdiff(temp_old$location_code, temp_new$location_code)
      # 
      # temp_old %>%
      #   dplyr::rename(location_old = location_code) %>%
      #   cbind(dplyr::select(temp_new, location_code)) %>%
      #   filter(location_old != location_code) %>%
      #   View()
      # 
      # 
      # temp_old %>%
      #   dplyr::rename(latitude_old = latitude,
      #          longitude_old = longitude) %>%
      #   cbind(dplyr::select(temp_new, latitude, longitude)) %>%
      #   filter(latitude != latitude_old | longitude != longitude_old) %>%
      #   View()

      # Location code 
      # 20722; 20875; 26866; 30325; 56148; 56151; 74393; 77384
      
    }
  }
}

# Make time trend plots---------------------------------------------------------------------------------

geomatch_time_plots <- function(geomatch_data,
                                old_data,
                                save_plot = paste0('/lbd_hiv/', "/data/covariate_surveys/", topic)) {
  # Color palette
  carto_discrete <- rep(c("#7F3C8D","#11A579","#F2B701","#E73F74",
                          "#3969AC","#80BA5A","#E68310","#008695",
                          "#CF1C90","#f97b72","#4b4b8f","#A5AA99"), 3)
  
  # Make data plots per each country
  gg_plots <- 
    lapply(sort(unique(geomatch_data$country)), function(cntry) {
      
      message(paste("Making time trend plot for", cntry))
      data <- geomatch_data[country == cntry]
      
      data_collapse <- 
        data[, .(indicator = weighted.mean(indicator, pweight, na.rm = T),
                 s_dev = sd(indicator, na.rm = T),
                 N_obs = .N,
                 N = sum(pweight) ^ 2 / sum(pweight ^ 2)), 
             by = .(nid, country, survey_name, year, point)][order(year)]
      
      val_range <- c(0, 1)
      
      gg_admin <- 
        ggplot(data = data_collapse, aes(x = year, y = indicator)) +
        geom_point(aes(size = N, shape = factor(point), fill = as.factor(str_trunc(survey_name, 10)))) +
        geom_smooth(se = F, color = "grey", alpha = 0.7, size = 0.4, linetype = "dashed", method = "lm") +
        theme_bw(base_size = 16) +
        scale_fill_manual("Survey", values = carto_discrete, drop = F) +
        scale_shape_manual("Type", 
                           values = c("0" = 22, "1" = 21, "2" = 12, "3" = 10),
                           label = c("0" = "polygon", "1" = "point", "2" = "New survey", "3" = "New survey"), drop = F) +
        coord_cartesian(ylim = val_range) +
        scale_x_continuous(limits = c(2000, 2018)) +
        scale_size(range = c(1,7)) +
        theme(
          strip.background = element_blank(),
          plot.caption = element_text(hjust = 0.5),
          plot.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 7),
          legend.justification = "top"
        ) +
        guides(fill = guide_legend(order = 1, override.aes = list(size = 5, shape = 21)),
               shape = guide_legend(order = 2, override.aes = list(size = 5))) +
        labs(y = paste(indicator, "prevalence"), x = "Year", title = cntry)
      return(gg_admin)
    })
  
  names(gg_plots) <- sort(unique(geomatch_data$country))
  
  # Make a pdf of plot of all countries
  message(paste("Saving plots of all countries to", paste0(save_plot, "/geomatch_trend_plots_all.pdf")))
  pdf(paste0(save_plot, "/geomatch_trend_plots_all_.pdf"), height = 11, width = 18)
  lapply(gg_plots, print)
  dev.off()
  
  # Now only plot countries that have added new data -----------------------------------------------------------------------------
  new_nid <- setdiff(geomatch_data$nid, old_data$nid)
  new_country <- 
    geomatch_data %>% 
    filter(nid %in% new_nid) %>% 
    pull(country) %>% unique() %>% sort()
  
  new_plots <- 
    lapply(new_country, function(cntry){
      
      message(paste("Making time trend plot for", cntry))
      data <- geomatch_data[country == cntry][, year := as.numeric(year)]
      
      # Right now drop pweight and point 
      data <- data[!is.na(point) & !is.na(pweight)]
      
      data_collapse <- 
        data[, .(indicator = weighted.mean(indicator, pweight, na.rm = T),
                 s_dev = sd(indicator, na.rm = T),
                 N_obs = .N,
                 N = sum(pweight) ^ 2 / sum(pweight ^ 2)), 
             by = .(nid, country, survey_name, year, point)][order(year)] %>% 
        mutate(point = as.numeric(point)) %>% 
        mutate(point = ifelse(nid %in% new_nid, point + 2, point))
      
      val_range <- c(0, 1)
      
      
      gg_admin <- 
        ggplot(data = data_collapse, aes(x = year, y = indicator)) +
        geom_point(aes(size = N, shape = as.factor(point), fill = as.factor(str_trunc(survey_name, 10)))) +
        geom_smooth(se = F, color = "grey", alpha = 0.7, size = 0.4, linetype = "dashed", method = "lm") +
        theme_bw(base_size = 16) +
        scale_shape_manual("Type", 
                           values = c("0" = 22, "1" = 21, "2" = 12, "3" = 10),
                           label = c("0" = "polygon", "1" = "point", "2" = "New survey", "3" = "New survey"), drop = F) +
        scale_fill_manual("Survey", values = carto_discrete, drop = F) +
        scale_color_manual("Survey", values = carto_discrete, drop = F) +
        coord_cartesian(ylim = val_range) +
        scale_x_continuous(limits = c(2000, 2018)) +
        scale_size(range = c(1,7)) +
        theme(
          strip.background = element_blank(),
          plot.caption = element_text(hjust = 0.5),
          plot.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 7),
          legend.justification = "top"
        ) +
        guides(fill = guide_legend(order = 1, override.aes = list(size = 5, shape = 21)),
               shape = guide_legend(order = 2, override.aes = list(size = 5))) +
        labs(y = paste(indicator, "prevalence"), x = "Year", title = cntry)
      return(gg_admin)
    })
  
  message(paste("Saving plots of new countries to", paste0(save_plot, "/geomatch_trend_plots_all.pdf")))
  pdf(paste0(save_plot, "/geomatch_trend_plots_new.pdf"), height = 11, width = 18)
  lapply(new_plots, print)
  dev.off()
  message("Done!")
  
}
