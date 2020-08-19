detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices",
                      "package:utils","package:datasets","package:methods",
                      "package:base", "package:rio")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) 
    detach(package, character.only=TRUE)
}
detachAllPackages()

library(dplyr)
library(data.table)

rm(list = ls())
indi_fam <- 'water'

#for some reason this column is a factor?
read_convert_weights <- function(x){
  temp <- fread(x, stringsAsFactors = F)
  temp$sum_of_sample_weights <- as.numeric(temp$sum_of_sample_weights)
  return(temp)
}

exc_sani <- fread('<<<< FILEPATH REDACTED >>>>')
exc_water <- fread('<<<< FILEPATH REDACTED >>>>')

format_id <- function(df, indicator) {
  df <- as.data.frame(df)
  remove_col <- grep('V', names(df))
  df <- df[,setdiff(names(df), names(df)[remove_col])]
  df <- df %>%
    group_by(nid) %>% 
    mutate(sum_of_sample_weights = ifelse(is.na(sum_of_sample_weights), 
                                          weight * N, 
                                          sum_of_sample_weights)) %>%
    mutate(year = weighted.mean(x = year, w = sum_of_sample_weights)) %>%
    ungroup() 
  df$N <- round(df$N, digits = 7)
  df[,indicator] <- round(df[,indicator], digits = 7)
  df <- filter(df, N > 0)
  df <- filter(df, year > 1999)
  df <- filter(df, !is.na(latitude), !is.na(longitude))
  df$year <- round(df$year)
  
  # recalculate prop
  df$prop <- as.numeric(df[[indicator]])/as.numeric(df$N)
  
  df <- subset(df, nid %in% exc_nids)
  df <- as.data.table(df)
  return(df)
}

if (indi_fam == 'water' ) {
  exc_nids <- exc_water$nid
  
  setwd('<<<< FILEPATH REDACTED >>>>')
  network <- lapply(list.files(), read_convert_weights)
  network <- do.call(rbind, network)
  
  setwd('<<<< FILEPATH REDACTED >>>>')
  imp <- lapply(list.files(), read_convert_weights)
  imp <- do.call(rbind, imp)
  
  setwd('<<<< FILEPATH REDACTED >>>>')
  piped <- lapply(list.files(), read_convert_weights)
  piped <- do.call(rbind, piped)
  
  setwd('<<<< FILEPATH REDACTED >>>>')
  unimp <- lapply(list.files(), read_convert_weights)
  unimp <- do.call(rbind, unimp)
  
  for (i in c('network','piped','imp', 'unimp')) {  
    mydat <- get(i)
    
    names(mydat)[which(names(mydat) == i)] <- 'indi'
    
    mydat <- mydat %>% dplyr::select(-lat.y, -long.y) %>%
      rename(latitude = lat.x, longitude = long.x,
             year = year_start,
             country = iso3) %>% 
      arrange(nid, country, year, shapefile, location_code) %>%
      mutate(indi_bin = indi*N, merge_id = 1:n()) %>%
      rename(indi_prop = indi)
    
    names(mydat)[which(names(mydat) == 'indi_bin')] <- paste0('w_',i)
    names(mydat)[which(names(mydat) == 'indi_prop')] <- paste0(i, '_prop')
    assign(i,mydat)
  }
  
  imp_denom <- dplyr::select(imp, w_imp, merge_id)
  imp_denom <- distinct(imp_denom)
  piped_denom <- dplyr::select(piped, w_piped, merge_id)
  piped_denom <- distinct(piped_denom)
  
  unimp <- left_join(unimp, imp_denom, by = 'merge_id')
  unimp <- left_join(unimp, piped_denom, by = 'merge_id')
  
  network <- left_join(network, piped_denom, by = 'merge_id')
  network <- mutate(network, N = w_piped) %>%
    rename(w_network_cr = network_prop, prop = network_prop) %>%
    mutate(w_network_cr = prop * N) %>%
    select(-w_piped, -w_network, -merge_id) %>%
    filter(N > 0)
  
  unimp <- mutate(unimp, N = ((N)) - (w_imp) - w_piped) %>% 
    mutate(unimp_prop = w_unimp/N) %>%
    rename(prop = unimp_prop, w_unimp_cr = w_unimp) %>%
    dplyr::select(-w_imp, -w_piped, -merge_id) %>%
    filter(N > 0)
  
  piped <- left_join(piped, imp_denom, by = 'merge_id')
  piped <- mutate(piped, piped_prop = w_piped/N) %>%
    rename(prop = piped_prop) %>%
    dplyr::select(-w_imp, -merge_id) %>%
    filter(N > 0)
  
  imp <- left_join(imp, piped_denom, by = 'merge_id')
  
  imp <- mutate(imp, N = N - w_piped) %>% 
    mutate(imp_prop = w_imp/N) %>%
    rename(prop = imp_prop, w_imp_cr = w_imp) %>%
    dplyr::select(-w_piped, -merge_id) %>%
    filter(N > 0)
  rm(mydat, imp_denom, piped_denom)
  
  piped_pt <- read.csv('<<<< FILEPATH REDACTED >>>>',
                       stringsAsFactors = F)
  piped_pt <- dplyr::select(piped_pt, -imp)
  piped <- rbind(piped, piped_pt)
  piped <- as.data.table(piped)
  piped <- format_id(piped, 'w_piped')
  write.csv(piped, '<<<< FILEPATH REDACTED >>>>', row.names = FALSE)
  rm(piped_pt)
  
  network_pt <- read.csv('<<<< FILEPATH REDACTED >>>>',
                         stringsAsFactors = F)

  network <- rbind(network, network_pt)
  network <- as.data.table(network)
  network[is.na(year_median), year_median := surv_year]
  network <- format_id(network, 'w_network_cr')
  write.csv(network, '<<<< FILEPATH REDACTED >>>>', row.names = FALSE)
  rm(network_pt)
  
  unimp_pt <- read.csv('<<<< FILEPATH REDACTED >>>>',
                       stringsAsFactors = F)

  unimp_pt <- select(unimp_pt, -piped,-w_piped)
  unimp <- rbind(unimp, unimp_pt)
  unimp <- format_id(unimp, 'w_unimp_cr')
  write.csv(unimp, '<<<< FILEPATH REDACTED >>>>', row.names = FALSE)
  rm(unimp_pt)
  
  imp_pt <- read.csv('<<<< FILEPATH REDACTED >>>>',
                     stringsAsFactors = F)

  imp <- rbind(imp, imp_pt)
  imp <- format_id(imp, 'w_imp_cr')
  write.csv(imp, '<<<< FILEPATH REDACTED >>>>', row.names = FALSE)
  rm(imp_pt)
}

indi_fam <- 'sani'
if (indi_fam == 'sani' ) {
  exc_nids <- exc_sani$nid
  
  setwd('<<<< FILEPATH REDACTED >>>>')
  network <- lapply(list.files(), read_convert_weights)
  network <- do.call(rbind, network)
  
  setwd('<<<< FILEPATH REDACTED >>>>')
  piped <- lapply(list.files(), read_convert_weights)
  piped <- do.call(rbind, piped)
  
  setwd('<<<< FILEPATH REDACTED >>>>')
  imp <- lapply(list.files(), read_convert_weights)
  imp <- do.call(rbind, imp)
  
  setwd('<<<< FILEPATH REDACTED >>>>')
  unimp <- lapply(list.files(),read_convert_weights)
  unimp <- do.call(rbind, unimp)
  
  for (i in c('network','piped', 'imp', 'unimp')) {  
    mydat <- get(i)
    
    names(mydat)[which(names(mydat) == i)] <- 'indi'
    
    mydat <- mydat %>% select(-lat.y, -long.y) %>%
      rename(latitude = lat.x, longitude = long.x,
             year = year_start,
             country = iso3) %>% 
      arrange(nid, country, year, shapefile, location_code) %>%
      mutate(indi_bin = indi*N, merge_id = 1:n()) %>%
      rename(indi_prop = indi)
    
    names(mydat)[which(names(mydat) == 'indi_bin')] <- paste0('s_',i)
    names(mydat)[which(names(mydat) == 'indi_prop')] <- paste0(i, '_prop')
    
    assign(i,mydat)
  }
  
  imp_denom <- select(imp, s_imp, merge_id)
  piped_denom <- select(piped, s_piped, merge_id)

  
  imp <- left_join(imp, piped_denom, by = 'merge_id')
  imp <- mutate(imp, N = N - s_piped) %>% 
    mutate(imp_prop = s_imp/N) %>%
    rename(s_imp_cr = s_imp, prop = imp_prop) %>%
    select(-s_piped, -merge_id) %>%
    filter(N > 0)
  
  unimp <- left_join(unimp, imp_denom, by = 'merge_id')
  unimp <- left_join(unimp, piped_denom, by = 'merge_id')
  unimp <- mutate(unimp, N = N - s_imp - s_piped) %>% 
    mutate(unimp_prop = s_unimp/N) %>%
    rename(s_unimp_cr = s_unimp, prop = unimp_prop) %>%
    select(-s_imp, -s_piped, -merge_id) %>%
    filter(N > 0)
  
  network <- left_join(network, piped_denom, by = 'merge_id')
  network <- mutate(network, N = s_piped) %>%
    rename(s_network_cr = network_prop, prop = network_prop) %>%
    mutate(s_network_cr = prop * N) %>%
    select(-s_piped, -s_network, -merge_id) %>%
    filter(N > 0)
  
  piped <- rename(piped, prop = piped_prop) %>% 
    select(-merge_id)
  rm(mydat, imp_denom, piped_denom)
  
  unimp_pt <- read.csv('<<<< FILEPATH REDACTED >>>>',
                       stringsAsFactors = F)
  unimp_pt <- select(unimp_pt, -piped, -network)
  unimp <- rbind(unimp, unimp_pt)
  unimp <- format_id(unimp, 's_unimp_cr')
  write.csv(unimp, '<<<< FILEPATH REDACTED >>>>', row.names = FALSE)
  rm(unimp_pt)
  
  network_pt <- read.csv('<<<< FILEPATH REDACTED >>>>',
                         stringsAsFactors = F)
  network <- rbind(network, network_pt)
  network <- format_id(network, 's_network_cr')
  write.csv(network, '<<<< FILEPATH REDACTED >>>>', row.names = FALSE)
  rm(network_pt)
  
  
  imp_pt <- read.csv('<<<< FILEPATH REDACTED >>>>',
                     stringsAsFactors = F)
  imp <- rbind(imp, imp_pt)
  imp <- format_id(imp, 's_imp_cr')
  write.csv(imp, '<<<< FILEPATH REDACTED >>>>', row.names = FALSE)
  rm(imp_pt)
  
  piped_pt <- read.csv('<<<< FILEPATH REDACTED >>>>',
                       stringsAsFactors = F)
  piped_pt <- select(piped_pt, -network)
  piped <- rbind(piped, piped_pt)
  piped <- format_id(piped, 's_piped')
  write.csv(piped, '<<<< FILEPATH REDACTED >>>>', row.names = FALSE)
  rm(piped_pt)
}
