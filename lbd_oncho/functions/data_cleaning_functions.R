### Data cleaning functions
###############################################################################################################

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
  
  # manual_check_summary <- unique(d[easy == 0, .(group_id, latitude, longitude, year, ihme_loc_id)])
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
  # return(list(obvious_drops, NULL, NULL))
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

split_dups_obvi_only_sample_size_off_by_one <- function(data, group_id, uid = "Master_UID", check_cols, calc_var_cols, get_unique, preferred_source = NULL, indicator = NULL) {
  d <- copy(data) %>% as.data.table()
  d[, group_id := get(group_id)]
  d[, m_uid := get(uid)]
  
  # vector of uids from duplicate clusters that should be dropped. Only 1 observation from each duplicate cluster is exempt from this and would be kept.
  obvious_drops <- c()
  
  # vector of uids from duplicate clusters that cannot be easily proessed and require more manual checks.
  manual_check <- c()
  
  for (i in unique(d$group_id)) {
    print(paste0(which(unique(d$group_id) == i), " of ", length(unique(d$group_id))))
    if ((nrow(unique(d[group_id == i, check_cols, with = F])) == 1) | (((max(unique(d[group_id == i, N])) - min(unique(d[group_id == i, N]))) %in% c(0, 1)) & (max(unique(d[group_id == i, ..indicator])) - min(unique(d[group_id == i, ..indicator]))) %in% c(0, 1))) {
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

plot_points <- function(data, lat_col = "latitude", long_col = "longitude", boundaries = "default", b_map, o_master, drop_outside_geo = F, fill_col = NULL, size_col = NULL, title_text) {
  if (boundaries == "default") {
    if (!("background_map" %in% ls())) {
      master_shape_all <- readRDS(<<<< FILEPATH REDACTED >>>>)
      gaul_list <- get_gaul_codes("lf_endem_afr")
      background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
    }
    
    if (!("outline_master" %in% ls())) {
      outline_master <- fast_load_shapefile("lf_g2015_2014_1")
    }
  } else {
    background_map <- b_map
    outline_master <- o_master
  }
  
  
  # A blank theme
  theme_empty <- theme_classic() +
    theme(
      axis.line = element_blank(), axis.text.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  
  background_extent <- extent(background_map)
  background_outline <- suppressMessages(
    geom_path(
      data = background_map,
      aes(
        x = long,
        y = lat,
        group = group
      ),
      color = "black",
      size = 1.4
    )
  )
  
  # adding admin boundaries to plots
  outline <- raster::crop(outline_master, background_extent)
  adm1_outline <- suppressMessages(
    geom_path(
      data = outline,
      aes(
        x = long,
        y = lat,
        group = group
      ),
      color = "black",
      linetype = 2
    )
  )
  
  background_map_not_endemic <- NULL
  background_map_overlay <- background_map
  background_map <- suppressMessages(
    geom_polygon(
      data = background_map,
      aes(
        x = long,
        y = lat,
        group = group
      ),
      fill = "white"
    )
  )
  
  d <- copy(data) %>% as.data.table()
  d[, latitude := get(lat_col)]
  d[, longitude := get(long_col)]
  
  if (drop_outside_geo) {
    ## Drop points/polygon out of background map extent
    d <- d[!(latitude < background_extent@ymin | latitude > background_extent@ymax |
               longitude < background_extent@xmin | longitude > background_extent@xmax), ]
  }
  
  ## Initialize map
  g_map <- ggplot() + background_map + background_map_not_endemic
  
  if (is.null(fill_col) == T & is.null(size_col) == T) {
    g_map <- g_map +
      adm1_outline +
      background_outline +
      geom_point(
        data = d,
        aes(
          x = longitude,
          y = latitude
        ),
        alpha = 1,
        size = 1.5,
        color = "black",
        fill = "red",
        shape = 21,
        stroke = 0.1
      ) +
      coord_equal() +
      theme_empty +
      ggtitle(title_text)
  } else {
    if (is.null(fill_col) == T) {
      g_map <- g_map +
        adm1_outline +
        background_outline +
        geom_point(
          data = d,
          aes(
            x = longitude,
            y = latitude,
            size = get(size_col)
          ),
          alpha = 1,
          color = "black",
          fill = "red",
          shape = 21,
          stroke = 0.1
        ) +
        coord_equal() +
        theme_empty +
        ggtitle(title_text) +
        scale_size_continuous(range = c(1.5, 4)) +
        guides(size = guide_legend(size_col))
    }
    if (is.null(size_col) == T) {
      g_map <- g_map +
        adm1_outline +
        background_outline +
        geom_point(
          data = d,
          aes(
            x = longitude,
            y = latitude,
            fill = get(fill_col)
          ),
          alpha = 1,
          size = 1.5,
          color = "black",
          shape = 21,
          stroke = 0.1
        ) +
        scale_fill_continuous(low = "yellow", high = "red") +
        coord_equal() +
        theme_empty +
        ggtitle(title_text) +
        guides(fill = guide_legend(fill_col))
    }
    if (is.null(fill_col) == F & is.null(size_col) == F) {
      g_map <- g_map +
        adm1_outline +
        background_outline +
        geom_point(
          data = d,
          aes(
            x = longitude,
            y = latitude,
            fill = get(fill_col),
            size = get(size_col)
          ),
          alpha = 1,
          color = "black",
          shape = 21,
          stroke = 0.1
        ) +
        scale_fill_continuous(low = "yellow", high = "red") +
        scale_size_continuous(range = c(1.5, 4)) +
        coord_equal() +
        theme_empty +
        ggtitle(title_text) +
        guides(size = guide_legend(size_col), fill = guide_legend(fill_col))
    }
  }
  
  
  return(g_map)
}

map_duplicates <- function(data, lat_col = "latitude", long_col = "longitude", boundaries = "default", drop_outside_geo = F, title) {
  d <- copy(data) %>% as.data.table()
  d[, latitude := get(lat_col)]
  d[, longitude := get(long_col)]
  
  locs <- unique(d[, .(latitude, longitude)])
  locs$loc_id <- c(1:nrow(locs))
  d <- merge(d, locs, by = c("latitude", "longitude"))
  if ("num_dups" %in% names(d)) {
    locs$num_records <- sapply(locs$loc_id, FUN = function(x) sum(d[loc_id == x, num_dups]))
  } else {
    locs$num_records <- sapply(locs$loc_id, FUN = function(x) nrow(d[loc_id == x, ]))
  }
  
  locs$num_years <- sapply(locs$loc_id, FUN = function(x) length(unique(d[loc_id == x, year])))
  
  map <- plot_points(data = locs, title_text = title, fill_col = "num_years", size_col = "num_records")
  return(map)
}

identify_wrong_adm0 <- function(data, s_poly = simple_polygon, fast_shapefiles = TRUE, rast = simple_raster) {
  library(stringr)
  d <- copy(data) %>% as.data.table()
  
  m_shape_cropped <- copy(s_poly)
  
  append <- d[is.na(latitude) | is.na(longitude)]
  
  d <- d[!is.na(latitude) & !is.na(longitude)]
  
  pts <- d[!is.na(latitude) & !is.na(longitude), .(latitude, longitude)]
  coordinates(pts) <- ~ longitude + latitude
  proj4string(pts) <- proj4string(m_shape_cropped)
  adm_data <- over(pts, m_shape_cropped) %>% as.data.table() %>% .[, .(ADM0_CODE)]
  d <- cbind(d, adm_data)
  
  ## Load link table for GADM codes
  adm_link <- fread(<<<< FILEPATH REDACTED >>>>)
  
  d <- merge(d, adm_link[, c("gadm_geoid", "iso3")], by.x = "ADM0_CODE", by.y = "gadm_geoid", all.x = TRUE)
  
  setnames(d, "iso3", "ihme_lc_id")
  
  d$country <- d$ihme_loc_id
  
  d$ihme_lc_id <- as.character(d$ihme_lc_id)
  d$country <- as.character(d$country)
  d[grep("_", d$ihme_lc_id), ihme_lc_id := unlist(str_split(ihme_lc_id, "_"))[[1]]]
  
  d$id <- 1:nrow(d)
  
  d_correct_country <- d[country == ihme_lc_id]
  
  d_wrong_country <- d[(country != ihme_lc_id) | is.na(ihme_lc_id)]
  
  cleaned_d <- copy(d_wrong_country)
  
  if (nrow(d_wrong_country) > 0) {
    message(paste0(nrow(d_wrong_country), " rows of data determined to be outside of their specified adm0"))
    message("Checking to see if any of these can be pushed to the appropriate adm0")
    
    updated_coords <- data.table()
    
    for (c in unique(d_wrong_country$country)) {
      print(c)
      country_shape <- m_shape_cropped[m_shape_cropped$ADM0_CODE == adm_link[iso3 == c, gadm_geoid],]
      if (length(country_shape) == 0) {
        print("Country not present in shapefile; skipping")
        next
      }
      country_rast <- rasterize(country_shape, rast)
      outside_pts <- d_wrong_country[country == c, .(longitude, latitude, id, ihme_lc_id)]
      land <- nearestLand(outside_pts[, 1:2], country_rast, 20000) # 20k buffer
      
      #check how many were moved
      print(sum(!is.na(land[, 1])))
      print(sum(is.na(land[, 1])))
      
      #create table with old and new coordinates for replacement
      pts_connect <- cbind(outside_pts, land)
      colnames(pts_connect) <- c("longitude_old", "latitude_old", "id", "ihme_lc_id", "longitude_new", "latitude_new")
      updated_coords <- rbind(updated_coords, pts_connect)
    }
    
    message(paste0(nrow(updated_coords[!is.na(latitude_new)]), " rows were nudged to the correct adm0."))
    message(paste0(nrow(updated_coords[is.na(latitude_new)]), " rows could not be nudged to the correct adm0. Please check these!"))
    
    # Make changes in dataset
    cleaned_d <- merge(cleaned_d, updated_coords, by = "id", all.x = TRUE)
    cleaned_d$latitude_new <- as.numeric(cleaned_d$latitude_new)
    cleaned_d$longitude_new <- as.numeric(cleaned_d$longitude_new)
    cleaned_d[, c("latitude", "longitude") := list(latitude_new, longitude_new)]
    
    # Deal with Sudan and South Sudan (special case)
    cleaned_d[ihme_loc_id == "SDN" & ihme_lc_id.x == "SSD" & !is.na(latitude_old) & !is.na(longitude_old), c("ihme_loc_id", "latitude", "longitude") := list("SSD", latitude_old, longitude_old)]
    
    d_wrong_country <- cleaned_d[is.na(latitude) | is.na(longitude)]
    cleaned_d <- cleaned_d[!is.na(latitude) & !is.na(longitude)]
  }
  
  append <- rbind(append, d_correct_country, fill = TRUE)
  return(list(rbind(append, cleaned_d, fill = TRUE), d_wrong_country))
}

create_check_dups <- function(check_country = NULL) {
  cyl <- unique(full_output[point == 1 & country == check_country, .(latitude_rounded, longitude_rounded, year, diagnostic)])
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

create_check_dups_no_year <- function(check_country = NULL) {
  cyl <- unique(full_output[point == 1 & country == check_country, .(latitude_rounded, longitude_rounded)])
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

create_check_dups_no_year_raw_coords <- function(check_country = NULL) {
  cyl <- unique(full_output[point == 1 & country == check_country, .(latitude, longitude, diagnostic)])
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

create_check_dups_no_year_raw_coords_location_name <- function(check_country = NULL) {
  cyl <- unique(full_output[point == 1 & country == check_country, .(latitude, longitude, location_name, diagnostic)])
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

create_check_dups_no_country <- function(d = full_output) {
  cyl <- unique(d[point == 1, .(latitude_rounded, longitude_rounded, year)])
  cyl$cyl.id <- c(1:nrow(cyl))
  temp_output <- as.data.table(left_join(d[point == 1], cyl, all.x = T))
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

create_check_dups_no_country_no_year <- function(d = full_output) {
  cyl <- unique(d[point == 1, .(latitude_rounded, longitude_rounded)])
  cyl$cyl.id <- c(1:nrow(cyl))
  temp_output <- as.data.table(left_join(d[point == 1], cyl, all.x = T))
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
