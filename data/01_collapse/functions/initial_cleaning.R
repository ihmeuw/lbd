initial_cleaning <- function(mydat = pt_collapse, var_family = indi_fam, 
                             dat_type = data_type, census = ipums) {
  
  message('Subset to relevant variables')
  if(!"int_year" %in% colnames(mydat)){
    mydat$int_year <- mydat$year_start
  }
  
  if (var_family == 'water') {
    ptdat_0 <- dplyr::select(mydat, nid, iso3, lat, long, survey_series, 
                             hhweight, urban, w_source_drink, hh_size, 
                             year_start, int_year, hhweight,
                             shapefile,location_code)
  } 

  if (var_family == 'sani') {
    if (census) {
      ptdat_0 <- dplyr::select(mydat, nid, iso3, lat, long, survey_series, 
                               hhweight, urban, t_type, shared_san, hh_size, 
                               year_start, int_year, hhweight, shapefile, 
                               location_code, sewage)  
    } else {
      ptdat_0 <- dplyr::select(mydat, nid, iso3, lat, long, survey_series, 
                               hhweight, urban, t_type, shared_san, hh_size, 
                               year_start, int_year,hhweight, shapefile, 
                               location_code)

    }
  }

  problem_list <- filter(ptdat_0, hh_size <= 0)
  
  ptdat_0$true_weight <- ptdat_0$hhweight
  
  message('Create a unique cluster id')
  if (dat_type == 'pt') {
    ptdat <- mutate(ptdat_0, 
                    cluster_id = paste(iso3, lat, long, nid, 
                                       year_start, sep = "_"))
  } else {
    ptdat <- mutate(ptdat_0, 
                    cluster_id = paste(iso3, shapefile, 
                                       location_code, nid, 
                                       year_start, sep = "_"))  
  }

  message('Create a table which assigns numbers to unique IDs 
          and merge it back to data to have shorter
          unique IDs')
  short_id <- data.frame(cluster_id = unique(ptdat$cluster_id), 
                         id_short = seq(1:length(unique(ptdat$cluster_id))),
                         stringsAsFactors = F)
  ptdat <- left_join(ptdat, short_id, by = 'cluster_id')
  rm(short_id)

  message('Remove longer cluster_ids')
  ptdat <- dplyr::select(ptdat, -cluster_id)

  message('Change weight to 1 if collapsing point data')
  if (dat_type == "pt" & agg_level != 'country') {ptdat$hhweight <- 1}

  message('Change shapefile and location code to missing if collapsing point data')
  if (dat_type == "pt") {ptdat$shapefile <- NA; ptdat$location_code <- NA}
  
  message('Removing points with not lats or longs')
  if (dat_type == "pt" & agg_level != 'country') {
    ptdat <- subset(ptdat, !is.na(lat) | !is.na(long))}
  
  message('Finding median interview year')
  ptdat <- ptdat %>% group_by(id_short) %>%
    mutate(int_year = median(int_year, na.rm = TRUE),
           year_median = weighted.mean(x = int_year, w = true_weight * hh_size)) %>% 
    ungroup

  results <- list(ptdat, ptdat_0)
  return(results)
}

