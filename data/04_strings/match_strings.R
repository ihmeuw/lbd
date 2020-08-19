# Clear environment
rm(list = ls())

package_lib <- '<<<< FILEPATH REDACTED >>>>'
.libPaths(package_lib)
library(magrittr)
library(data.table)
library(feather)

stg_mast <- fread('<<<< FILEPATH REDACTED >>>>')

most_recent_extracts <- list.files('<<<< FILEPATH REDACTED >>>>', 
                                   full.names = T, pattern=".feather$") %>% 
  grep(value=T, pattern="poly", invert=T)
extract_info <- file.info(most_recent_extracts)
extract_info$path <- rownames(extract_info)
extract_info <- data.table(extract_info)
most_recent_point <- extract_info[order(mtime, decreasing = T), path][1]

most_recent_extracts <- list.files('<<<< FILEPATH REDACTED >>>>', 
                                   full.names = T, pattern=".feather$") %>% 
  grep(value=T, pattern="point", invert=T)
extract_info <- file.info(most_recent_extracts)
extract_info$path <- rownames(extract_info)
extract_info <- data.table(extract_info)
most_recent_poly <- extract_info[order(mtime, decreasing = T), path][1]

message("Loading big extraction .Rdata")

point <- read_feather(most_recent_point)
pt1_collapse <- as.data.table(read_feather(most_recent_poly))
packaged <- as.data.table(rbind(point, pt_collapse))

packaged$w_source_drink <- tolower(packaged$w_source_drink)
packaged$w_source_other <- tolower(packaged$w_source_other)

packaged$t_type <- tolower(packaged$t_type)
packaged$sewage <- tolower(packaged$sewage)

packaged[!is.na(sewage) & !is.na(t_type), t_type := paste0(t_type, ' ', sewage)]
packaged[is.na(t_type) & !is.na(sewage), t_type := sewage]
packaged[,t_type := trimws(t_type)]
packaged[,w_source_drink := trimws(w_source_drink)]

packaged <- unique(packaged[,c('nid','iso3','year_start',
                               'w_source_drink','t_type')])
packaged <- subset(packaged, !iso3 %in% subset(stg_mast, Stage == 3)$iso3) %>% 
  subset(year_start > 1999)


w <- fread('<<<< FILEPATH REDACTED >>>>')
t <- fread('<<<< FILEPATH REDACTED >>>>')


w_new <- merge(unique(packaged[,c('nid','iso3','year_start','w_source_drink')]), w, 
               by = c('nid','iso3','year_start','w_source_drink'), 
               all.x = T, all.y = T)
t_new <- merge(unique(packaged[,c('nid','iso3','year_start','t_type')]), t, 
               by = c('nid','iso3','year_start','t_type'), all.x = T, all.y = T)


write.csv(w_new, '<<<< FILEPATH REDACTED >>>>', row.names = F)
write.csv(t_new, '<<<< FILEPATH REDACTED >>>>', row.names = F)
