### LF Data cleaning functions
###
#################################################################################################################


split_dups <- function(data, group_id, uid = "Master_UID", check_cols, calc_var_cols, get_unique, preferred_source = NULL) { # check_cols = c('N', 'had_lf_poly')
  d <- copy(data) %>% as.data.table()
  d[, group_id := get(group_id)]
  d[, m_uid := get(uid)]
  
  # vector of uids from duplicate clusters that should be dropped. Only 1 observation from each duplicate cluster is exempt from this and would be kept.
  obvious_drops <- c()
  
  # vector of uids from duplicate clusters that cannot be easily proessed and require more manual checks.
  manual_check <- c()
  
  for (i in unique(d$group_id)) {
    print(paste0(which(unique(d$group_id) == i), " of ", length(unique(d$group_id))))
    if (nrow(unique(d[group_id == i, check_cols, with = F])) == 1) {
      d[group_id == i, easy := 1]
      temp <- d[group_id == i,]
      if(!is.null(preferred_source)){
        if(preferred_source %in% unique(d[group_id == i, source])){ #if preferred source exists in the duplicate group, keep that data
          new_drops <- d[group_id == i & source != preferred_source, m_uid]
          if(nrow(d[group_id == i & source == preferred_source]) > 1){ #if there is duplication within preferred source, drop one
            new_drops <- c(new_drops, d[group_id == i & source == preferred_source, m_uid][2:length(d[group_id == i & source == preferred_source, m_uid])])
          }
          
        }else{#preferred source isn't available so no preference for which to drop
          new_drops <-  d[group_id == i, m_uid][2:length(d[group_id == i, m_uid])]
        }
      }else{#no preference for which one to drop/retain
        new_drops <-  d[group_id == i, m_uid][2:length(d[group_id == i, m_uid])]
      }
      
      obvious_drops <- c(obvious_drops, new_drops)
    } else {
      d[group_id == i, easy := 0]
      manual_check <- c(manual_check, d[group_id == i, m_uid])
    }
  }
  
  message(paste0(length(unique(d[easy == 1, group_id])), " duplicate clusters have the same values in check_cols and can be easily de-duplicated. ", length(unique(d[easy == 0, group_id])), " require manual checks."))
  
  manual_check_summary <- unique(d[easy == 0, .(group_id, latitude, longitude, year, country)])
  manual_check_summary$num_dups <- sapply(manual_check_summary$group_id, FUN = function(x) length(d[group_id == x, m_uid]))
  manual_check_summary$uids <- sapply(manual_check_summary$group_id, FUN = function(x) paste(d[group_id == x, m_uid], collapse = ";"))

  if (nrow(manual_check_summary) > 0) {
    calc_vars <- lapply(calc_var_cols, FUN = function(var_col) sapply(manual_check_summary$group_id, FUN = function(x) sqrt(var(d[group_id == x, var_col, with = F]))))
    for (i in c(1:length(calc_vars))) {
      manual_check_summary <- cbind(manual_check_summary, calc_vars[[i]])
      setnames(manual_check_summary, "V2", paste0("std.dev_", calc_var_cols[i]))
    }

    get_unique_vars <- lapply(get_unique, FUN = function(col) sapply(manual_check_summary$group_id, FUN = function(x) paste0(unlist(unique(d[group_id == x, col, with = F])), collapse = "; ")))
    if (length(get_unique_vars) > 0) {
      for (i in c(1:length(get_unique_vars))) {
        manual_check_summary <- cbind(manual_check_summary, get_unique_vars[[i]])
        setnames(manual_check_summary, "V2", paste0("unique_", get_unique[i]))
      }
    }
  }
  
  return(list(obvious_drops, manual_check, manual_check_summary))
}

split_dups_obvi_only <- function(data, group_id, uid = "Master_UID", check_cols, calc_var_cols, get_unique, preferred_source = NULL) { # check_cols = c('N', 'had_lf_poly')
  d <- copy(data) %>% as.data.table()
  d[, group_id := get(group_id)]
  d[, m_uid := get(uid)]
  
  # vector of uids from duplicate clusters that should be dropped. Only 1 observation from each duplicate cluster is exempt from this and would be kept.
  obvious_drops <- c()
  
  # vector of uids from duplicate clusters that cannot be easily proessed and require more manual checks.
  manual_check <- c()
  
  for (i in unique(d$group_id)) {
    print(paste0(which(unique(d$group_id) == i), " of ", length(unique(d$group_id))))
    if (nrow(unique(d[group_id == i, check_cols, with = F])) == 1) {
      d[group_id == i, easy := 1]
      temp <- d[group_id == i,]
      if(!is.null(preferred_source)){
        if(preferred_source %in% unique(d[group_id == i, source])){ #if preferred source exists in the duplicate group, keep that data
          new_drops <- d[group_id == i & source != preferred_source, m_uid]
          if(nrow(d[group_id == i & source == preferred_source]) > 1){ #if there is duplication within preferred source, drop one
            new_drops <- c(new_drops, d[group_id == i & source == preferred_source, m_uid][2:length(d[group_id == i & source == preferred_source, m_uid])])
          }
          
        }else{#preferred source isn't available so no preference for which to drop
          new_drops <-  d[group_id == i, m_uid][2:length(d[group_id == i, m_uid])]
        }
      }else{#no preference for which one to drop/retain
        new_drops <-  d[group_id == i, m_uid][2:length(d[group_id == i, m_uid])]
      }
      
      obvious_drops <- c(obvious_drops, new_drops)
    } else {
      d[group_id == i, easy := 0]
    }
  }
  
  message(paste0(length(unique(d[easy == 1, group_id])), " duplicate clusters have the same values in check_cols and can be easily de-duplicated. ", length(unique(d[easy == 0, group_id])), " require manual checks."))
  
  return(list(obvious_drops, NULL, NULL))
}

create_check_dups <- function(check_country = NULL) {
  cyl <- unique(full_output[point == 1 & country == check_country, .(latitude_rounded, longitude_rounded, year)])
  cyl$cyl.id <- c(1:nrow(cyl))
  temp_output <- as.data.table(left_join(full_output[point == 1], cyl, all.x = T))
  cyl.tab <- temp_output$cyl.id %>%
    table() %>%
    as.data.table()
  setnames(cyl.tab, ".", "cyl.id")
  check_dups <- cyl.tab[N > 1]
  check_dups <- temp_output[cyl.id %in% check_dups$cyl.id]
  check_dups[, count := uniqueN(source), by = cyl.id]
  check_dups <- check_dups[count > 1]
  return(check_dups)
}

create_check_dups_rounded_1 <- function(check_country = NULL) {
  cyl <- unique(full_output[point == 1 & country == check_country, .(latitude_rounded_1, longitude_rounded_1, year)])
  cyl$cyl.id <- c(1:nrow(cyl))
  temp_output <- as.data.table(left_join(full_output[point == 1], cyl, all.x = T))
  cyl.tab <- temp_output$cyl.id %>%
    table() %>%
    as.data.table()
  setnames(cyl.tab, ".", "cyl.id")
  check_dups <- cyl.tab[N > 1]
  check_dups <- temp_output[cyl.id %in% check_dups$cyl.id]
  check_dups[, count := uniqueN(source), by = cyl.id]
  check_dups <- check_dups[count > 1]
  return(check_dups)
}

create_check_dups_same_source <- function(check_country = NULL) {
  cyl <- unique(full_output[point == 1 & country == check_country, .(latitude_rounded, longitude_rounded, year)])
  cyl$cyl.id <- c(1:nrow(cyl))
  temp_output <- as.data.table(left_join(full_output[point == 1], cyl, all.x = T))
  cyl.tab <- temp_output$cyl.id %>%
    table() %>%
    as.data.table()
  setnames(cyl.tab, ".", "cyl.id")
  check_dups <- cyl.tab[N > 1]
  check_dups <- temp_output[cyl.id %in% check_dups$cyl.id]
  return(check_dups)
}

create_check_dups_same_source_no_year <- function(check_country = NULL) {
  cyl <- unique(full_output[point == 1 & country == check_country, .(latitude_rounded, longitude_rounded)])
  cyl$cyl.id <- c(1:nrow(cyl))
  temp_output <- as.data.table(left_join(full_output[point == 1], cyl, all.x = T))
  cyl.tab <- temp_output$cyl.id %>%
    table() %>%
    as.data.table()
  setnames(cyl.tab, ".", "cyl.id")
  check_dups <- cyl.tab[N > 1]
  check_dups <- temp_output[cyl.id %in% check_dups$cyl.id]
  return(check_dups)
}

create_check_dups_full_coords_no_year <- function(check_country = NULL) {
  cyl <- unique(full_output[point == 1 & country == check_country, .(latitude, longitude)])
  cyl$cyl.id <- c(1:nrow(cyl))
  temp_output <- as.data.table(left_join(full_output[point == 1], cyl, all.x = T))
  cyl.tab <- temp_output$cyl.id %>%
    table() %>%
    as.data.table()
  setnames(cyl.tab, ".", "cyl.id")
  check_dups <- cyl.tab[N > 1]
  check_dups <- temp_output[cyl.id %in% check_dups$cyl.id]
  check_dups[, count := uniqueN(source), by = cyl.id]
  check_dups <- check_dups[count > 1]
  return(check_dups)
}

create_check_dups_poly <- function(check_country = NULL) {
  cyl <- unique(full_output[point == 0 & country == check_country, .(year, shapefile, location_code)])
  cyl$cyl.id <- c(1:nrow(cyl))
  temp_output <- as.data.table(left_join(full_output[point == 0], cyl, all.x = T))
  cyl.tab <- temp_output$cyl.id %>%
    table() %>%
    as.data.table()
  setnames(cyl.tab, ".", "cyl.id")
  check_dups <- cyl.tab[N > 1]
  check_dups <- temp_output[cyl.id %in% check_dups$cyl.id]
  check_dups[, count := uniqueN(source), by = cyl.id]
  check_dups <- check_dups[count > 1]
  return(check_dups)
}
