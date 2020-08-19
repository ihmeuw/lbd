rm_miss <- function(mydat = ptdat, var_family = indi_fam, 
                    agg = agg_level, dt_type = data_type) {
  
  if (var_family == 'water') {
    mydat <- rename(mydat, indi = piped)
  }
  
  if (var_family == 'sani') {
    mydat <- rename(mydat, indi = od)
  }
  
  # Calculate data missingness by cluster
  missing <- mutate(mydat, miss = ifelse(is.na(mydat$indi), 1, 0))
  missing <- mutate(missing, miss_wt = miss*hh_size)
  missing <- missing %>% group_by(id_short) %>% 
    summarise(pct_miss = sum(miss_wt)/sum(hh_size))
  # Remove clusters with more than 20% weighted missingness
  miss_clusters <- dplyr::select(filter(missing, pct_miss > 0.2), id_short)
  mydat <- filter(mydat, !(id_short %in% miss_clusters$id_short))
  
  if(agg == 'country' & dt_type == 'pt') {
    # Calculate weight missingness by cluster for country agg & pts
    missing <- mutate(mydat, miss = ifelse(is.na(mydat$hhweight), 1, 0))
    missing <- mutate(missing, miss_wt = miss*hh_size)
    missing <- missing %>% group_by(id_short) %>% 
      summarise(pct_miss = sum(miss_wt)/sum(hh_size))
  
    miss_clusters <- select(filter(missing, pct_miss > 0.2), id_short)
    mydat <- filter(mydat, !(id_short %in% miss_clusters$id_short))
  }
  
  if (var_family == 'water') {
    mydat <- rename(mydat, piped = indi)
  }
  
  if (var_family == 'sani') {
    mydat <- rename(mydat, od = indi)
  }

  return(mydat)
}
  
# ----Missingness-------------------------------------------------------------
#function to identify missingness in key variables, with option for weighted %
idMissing <- function(input.dt, this.var, criteria=.2, 
                      wt.var=NA, check.threshold=F, 
                      threshold=NA, debug=F) {
  
  #allow for interactive debugs
  if (debug) browser()
  
  #set as copy so you dont save the new vars
  dt <- copy(input.dt)
  
  #set the original row count in order to print the data loss due to missingness
  og.count <- nrow(dt)
  
  #define the weight variable 
  if (is.na(wt.var)) {
    dt[, wt := 1]
  } else dt[, wt := get(wt.var)]
  
  # Calculate data missingness by cluster for given variable
  # Alternatively check for validity against a minimum threshold (usually 0)
  if (!check.threshold) {
    dt[, miss := lapply(.SD, is.na), .SDcols=this.var]
  } else dt[, miss := lapply(.SD, function(x) x <= threshold), .SDcols=this.var]
  
  #calc pct miss, weight by selected variable (typically hh_size)
  dt[, pct_miss := sum(miss*wt, na.rm=T)/sum(wt, na.rm=T), by=id_short] 
  
  # Return the IDs of any clusters with > specified criteria weighted missingnesss or invalidity
  message("\nidentified #", dt[pct_miss>criteria, id_short] %>% 
            uniqueN, " clusters...or ~", 
          round((nrow(dt[pct_miss>criteria])/og.count)*100), 
          "% of rows \n based on criteria of >", criteria*100, "% ",
          ifelse(!check.threshold, 'missingness', 'invalidity'),
          " of ", this.var)
  
  clusters <- dt[pct_miss>criteria, id_short] %>% 
    unique
  
  if (length(clusters>1)) {  
    #save a diagnostic file with the clusters and the type of missingness
    dt <- dt[id_short %in% clusters, .(nid, iso3, int_year, id_short)] %>%
      setkey(., nid, iso3, int_year, id_short) %>% 
      unique(., by=key(.))
    
    dt[, count := sum(id_short %in% clusters), by=nid]
    dt[, var := this.var]
    dt[, type := ifelse(!check.threshold, 'missingness', 'invalidity')]
    return(clusters)
    
  } else return(NULL)
  
}