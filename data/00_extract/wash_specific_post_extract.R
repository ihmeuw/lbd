all[is.na(hhweight) & !is.na(pweight), hhweight := pweight]

drop <- c("line_id", "sex_id", "age_year", "age_month", "age_day", 
          "pweight", "cooking_fuel_mapped", "w_treat", "w_filter", 
          "w_boil", "w_bleach", "nid_n", "year", "t_type_multi", 
          "w_solar", "w_cloth", "w_settle", "hw_soap1", "hw_soap2", 
          "hw_soap3")
all <- all[, (drop):= NULL]

message("drop duplicate WN entries")
wn <- all[survey_module == "WN", ]
wn_key <- c("psu", "hh_id")
wn <- distinct(wn, psu, hh_id, .keep_all=T)
all <- all[survey_module != "WN", ]
all <- rbind(all, wn, fill=T, use.names=T)

message("custom wash fixes")
drop <- c("latitude", "longitude")
all <- all[, (drop):=NULL]

message("drop duplicate HH entries and cleanup hh_sizes")
#subset cases where all hh_sizes are present. Make sure each Row is a HH
has_hh_size_no_id <- all[!is.na(hh_size) & is.na(hh_id), ]
has_hh_size_id <- all[!is.na(hh_size) & !is.na(hh_id), ]
has_hh_size_id[, uq_id := paste(nid, psu, geospatial_id, 
                                hh_id, year_start, lat, long, 
                                shapefile, location_code, 
                                sep="_")] #includes space-time
has_hh_size_id[, prev_uq_id := paste(nid, psu, hh_id, sep="_")]
diff <- length(unique(has_hh_size_id$uq_id)) - length(unique(has_hh_size_id$prev_uq_id))
message(paste("There are", diff, 
              "more unique households from including spacetime than excluding."))
hhhs <- distinct(has_hh_size_id, uq_id, .keep_all=T)

#subset cases where all hh_sizes are missing and each row is not a HH. 
#   Set hh_size to 1
missing_hh_size <- all[is.na(hh_size) & survey_module != 'HH', ]
missing_hh_size[, hh_size := 1]

missing_hh_size_hh <- all[is.na(hh_size) & survey_module == 'HH', ]

packaged <- rbind(hhhs, has_hh_size_no_id, fill=T)
rm(hhhs)
rm(has_hh_size_no_id)
packaged <- rbind(packaged, missing_hh_size, fill=T)
rm(missing_hh_size)
packaged <- rbind(packaged, missing_hh_size_hh, fill=T)
rm(missing_hh_size_hh)

pt_collapse <- packaged[!is.na(lat) & !is.na(long), ]
poly_collapse <- packaged[(is.na(lat) | is.na(long)) & !is.na(shapefile) 
                          & !is.na(location_code), ]

library(feather)
message("Point Feather")
write_feather(pt_collapse, '<<<< FILEPATH REDACTED >>>>')
message("Poly Feather")
write_feather(poly_collapse, '<<<< FILEPATH REDACTED >>>>')