.libPaths('<<<< FILEPATH REDACTED >>>>')
rm(list = ls())

library(feather)
library(dplyr)

setwd('<<<< FILEPATH REDACTED >>>>')

indi_fam <- 'water'
if (indi_fam == 'water') {
  ptdat <- read_feather('<<<< FILEPATH REDACTED >>>>')
  ptdat <- filter(ptdat, !(is.na(lat)))
  
  w_network_cr <- ptdat
  w_network_cr <- dplyr::select(w_network_cr, -unimp, -surface)
  w_network_cr <- mutate(w_network_cr, point = 1, weight = 1, 
                         w_network_cr = (network*N*piped), N = N*piped)
  w_network_cr <- rename(w_network_cr, country = iso3, year = year_start, 
                         latitude = lat, longitude = long) %>% 
    dplyr::select(-imp,-piped, -network)
  w_network_cr <- filter(w_network_cr, N > 0)
  w_network_cr <- mutate(w_network_cr, prop = w_network_cr / N)
  write.csv(w_network_cr, file = '<<<< FILEPATH REDACTED >>>>', 
            row.names = FALSE)
  
  w_piped <- ptdat
  w_piped <- dplyr::select(w_piped, -unimp, -surface, -network)
  w_piped <- mutate(w_piped, point = 1, weight = 1, w_piped = (piped*N))
  w_piped <- rename(w_piped, country = iso3, year = year_start, prop = piped, 
                    N = N, latitude = lat,
                    longitude = long) %>%
    filter(N > 0)
  write.csv(w_piped, file = '<<<< FILEPATH REDACTED >>>>', 
            row.names = FALSE)
  
  w_imp_cr <- ptdat
  w_imp_cr <- dplyr::select(w_imp_cr, -surface, -unimp, -network)
  w_imp_cr <- mutate(w_imp_cr, point = 1, weight = 1, w_imp_cr = (imp*N), 
                     w_piped = (piped*N))
  w_imp_cr <- rename(w_imp_cr, country = iso3, year = year_start, N = N, 
                     latitude = lat,longitude = long)
  w_imp_cr <- mutate(w_imp_cr, N = ((N)) - (w_piped)) %>% 
    mutate(prop = w_imp_cr/N) %>%
    dplyr::select(-piped, -w_piped, -imp) %>%
    filter(N > 0)
  write.csv(w_imp_cr, '<<<< FILEPATH REDACTED >>>>', 
            row.names = FALSE)
  
  w_unimp_cr <- ptdat
  w_unimp_cr <- dplyr::select(w_unimp_cr, -surface, -network)
  w_unimp_cr <- mutate(w_unimp_cr, point = 1, weight = 1, w_unimp_cr = (unimp*N), 
                       w_imp = (imp*N), w_piped = (piped*N))
  w_unimp_cr <- rename(w_unimp_cr, country = iso3, year = year_start, N = N, 
                       latitude = lat,longitude = long)
  w_unimp_cr <- mutate(w_unimp_cr, N = ((N)) - (w_imp) - w_piped) %>% 
    mutate(prop = w_unimp_cr/N) %>%
    dplyr::select(-imp, -w_imp, -unimp, w_piped, piped) %>%
    filter(N > 0)
  write.csv(w_unimp_cr, '<<<< FILEPATH REDACTED >>>>', 
            row.names = FALSE)
  
}

rm(list = ls())
indi_fam <- 'sani'
if (indi_fam == 'sani') {
  ptdat <- read_feather('<<<< FILEPATH REDACTED >>>>')
  ptdat <- filter(ptdat, !(is.na(lat)))
  
  s_network_cr <- ptdat
  s_network_cr <- dplyr::select(s_network_cr, -unimp, -od)
  s_network_cr <- mutate(s_network_cr, point = 1, weight = 1, 
                         s_network_cr = (network*N*piped), N = N*piped)
  s_network_cr <- rename(s_network_cr, country = iso3, year = year_start, 
                         latitude = lat, longitude = long) %>% 
    dplyr::select(-imp,-piped,-network)
  s_network_cr <- filter(s_network_cr, N > 0)
  s_network_cr <- mutate(s_network_cr, prop = s_network_cr / N)
  write.csv(s_network_cr, file = '<<<< FILEPATH REDACTED >>>>', 
            row.names = FALSE)
  
  s_piped <- ptdat
  s_piped <- dplyr::select(s_piped, -od, -unimp,-imp)
  s_piped <- mutate(s_piped, point = 1, weight = 1, s_piped = (piped*N))
  s_piped <- rename(s_piped, country = iso3, year = year_start, prop = piped, 
                    N = N, latitude = lat, longitude = long)
  s_piped <- mutate(s_piped, N = (N))
  write.csv(s_piped, '<<<< FILEPATH REDACTED >>>>', 
            row.names = FALSE)
  
  s_imp_cr <- ptdat
  s_imp_cr <- dplyr::select(s_imp_cr, -od, -unimp, -network)
  s_imp_cr <- mutate(s_imp_cr, point = 1, weight = 1, s_imp_cr = (imp*N), 
                     s_piped = (piped*N))
  s_imp_cr <- rename(s_imp_cr, country = iso3, year = year_start, N = N, 
                     latitude = lat, longitude = long)
  s_imp_cr <- mutate(s_imp_cr, N = ((N)) - (s_piped)) %>% 
    mutate(prop = s_imp_cr/N) %>%
    dplyr::select(-piped, -s_piped, -imp) %>%
    filter(N > 0)
  write.csv(s_imp_cr, '<<<< FILEPATH REDACTED >>>>', 
            row.names = FALSE)
  
  s_unimp_cr <- ptdat
  s_unimp_cr <- dplyr::select(s_unimp_cr, -od)
  s_unimp_cr <- mutate(s_unimp_cr, point = 1, weight = 1, s_unimp_cr = (unimp*N), 
                       s_imp = (imp*N))
  s_unimp_cr <- rename(s_unimp_cr, country = iso3, year = year_start, N = N, 
                       latitude = lat, longitude = long)
  s_unimp_cr <- mutate(s_unimp_cr, N = ((N) - s_imp)) %>% 
    mutate(prop = s_unimp_cr/N) %>% 
    dplyr::select(-s_imp, -imp, -unimp) %>% 
    filter(N > 0)
  write.csv(s_unimp_cr, '<<<< FILEPATH REDACTED >>>>', 
            row.names = FALSE)
  
}