# Function to prep mortality data to be used in HIV small area estimation models 

#' @title vr_prep_hiv_mortality
#' @description: This function preps VR mortality data to be used in HIV small area estimation models. 
#' 
#'
#' @param country iso3 code for the country that you wish to prep HIV mortality data on
#' @param admin Admin level (for example "adm2")
#' @param years Years of population estimates 
#' @param shapefile_path Path to shapefile used 
#' 
#' @return Returns a data table that has the each admin unit by year and the number of HIV deaths after redistribution

vr_prep_hiv_mortality <- function(country, 
                                  admin,
                                  years) {
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(assertthat)
  
  # Link between GBD age group and SAE ages
  age_link <- 
    get_age_metadata(age_group_set_id = 12, gbd_round_id = 6) %>% 
    dplyr::select(age_group_id, age_group_years_start) %>% 
    dplyr::rename(age = age_group_years_start)
  
  # Get hiv cause id's
  # Maybe ammend these IDS
  source("<<<< FILEPATH REDACTED >>>>/get_cause_metadata.R")
  hiv_cause_ids <-
    get_cause_metadata(cause_set_id = 3, gbd_round_id = 6) %>%
    filter(acause_parent == "hiv") %>% pull(path_to_top_parent) %>% 
    str_split(",") %>% 
    unlist() %>% 
    unique() %>% 
    as.numeric() %>% 
    setdiff(c(294, 295, 955)) # Remote these higher up ones for now
  
  # Link between admin code and GBD location id
  core_repo <- paste0("<<<< FILEPATH REDACTED >>>>/lbd_core/")
  source(paste0(core_repo, "data_central/vr_functions.R"))
  
  # Link between admin2 code and GBD location ID 
  area_to_loc_id <-
    vr_pull_loc_meta(paste0(country, "_", admin)) %>%
    filter(level == as.numeric(str_extract(admin, "[0-9]+"))) %>%
    dplyr::select(area = uid, location_id) %>%
    unique()
  
  # Load country specific raw mortlity data (this may take some time)
  mort_raw <- vr_pull_cod(paste0(country,"_", admin), "redistribution")
  
  # Subset to HIV deaths
  mort_raw <- mort_raw[cause_id %in% hiv_cause_ids,]
  
  # Subset to correct years
  mort_raw <- mort_raw[year_id %in% years,]
  
  # See what percent of HIV deaths are mapped to admin1 instead of admin2
  admin1_code <- 
    vr_pull_loc_meta(paste0(country, "_", admin)) %>% 
    filter(level == 1) %>% 
    pull(location_id)
  
  pct_admin1 <-
    mort_raw %>%
    mutate(level = ifelse(location_id %in% admin1_code, 1, 2)) %>%
    group_by(level) %>%
    summarize(deaths = sum(deaths)) %>%
    ungroup() %>%
    mutate(pct = round(100*deaths / sum(deaths), 1)) %>%
    data.table()
  
  message(paste0(
      if (nrow(filter(mort_raw, location_id %in% admin1_code)) != 0) {
        pct_admin1[level == 1, pct]
      } else {
        0
      }, "% of deaths were mapped to admin 1 instead of admin 2.\nCurrently those are being dropped"))
  
  
  # Get correct age grouping and area ID for HIV  ----------------------------------------
  # For now omit things mapped to admin 1
  # Take raw mortality data and collapse by year and admin unit, keep age in 5 year bins from 0 to 80+
  mort_data <- 
    mort_raw %>%
    left_join(age_link, by = "age_group_id") %>% # get SAE age groupings, five years up until 80+
    mutate(age = case_when(age < 5 ~ 0,
                           between(age, 5, 75) ~ age,
                           age >= 80 ~ 80)) %>% # Get correct age mapping
    left_join(area_to_loc_id, by = "location_id") %>% # change to area
    dplyr::select(area, year_id, sex_id, age, deaths) %>% 
    dplyr::rename(year = year_id, sex = sex_id, events = deaths) %>% 
    group_by(area, year, sex, age) %>% # collapse by uid, year, sex, age
    dplyr::summarize(events = sum(events)) %>% 
    ungroup() %>% 
    filter(!is.na(area)) %>% # Temporary for now dropping deaths mapped to admin 1
    arrange(area, year, sex, age, events) %>% 
    data.table()
  
  
  # Add on admin-1 if required 
  if (nrow(mort_raw[location_id %in% admin1_code]) != 0) { 
    message("Assigning extra deaths to admin-1 to admin-2 proportional to deaths at that level")
    admin1_area_map <- 
      vr_pull_loc_meta(paste0(country, "_", admin)) %>% 
      filter(level == 2, year_end == max(year_end)) %>% # Gets most recent map to 
      dplyr::select(area = uid, admin2 = location_id, admin1 = location_parent_id) %>% 
      data.table()
    
    # Grab all deaths from raw data that are coded to admin-1
    mort_data_ad1 <-
      mort_raw %>% 
      filter(location_id %in% admin1_code) %>% 
      left_join(age_link, by = c("age_group_id")) %>% 
      mutate(age = case_when(age < 5 ~ 0,
                             between(age, 5, 75) ~ age,
                             age >= 80 ~ 80)) %>% 
      dplyr::select(admin1 = location_id, year = year_id, sex = sex_id, age, deaths) %>% 
      group_by(admin1, year, sex, age) %>% 
      dplyr::summarize(deaths = sum(deaths)) %>% 
      data.table()
    
    # Redistribute deaths for each area-year-age-sex proportional to how many deaths are already coded there
    mort_data_ad1_area <-
      rbindlist(mclapply(split(mort_data_ad1, 1:nrow(mort_data_ad1)), mc.cores = 25, function(x) {
        area_specific <-
          mort_data[area %in% admin1_area_map[admin1 == x[["admin1"]], area] & year == x[["year"]] & sex == x[["sex"]] & age == x[["age"]]]
        if (nrow(area_specific[events != 0]) != 0) {
          area_specific <-
            area_specific %>%
            mutate(pct = events / sum(events)) %>%
            mutate(deaths = x[["deaths"]]) %>%
            mutate(events = deaths * pct) %>%
            dplyr::select(area, year, sex, age, events) %>% 
            data.table()
        } else {
          # If there were no deaths recorded in area-year-age-sex then split equally
          area_specific <-
            area_specific %>%
            mutate(events = x[["deaths"]] / nrow(area_specific)) %>% 
            dplyr::select(area, year, sex, age, events) %>% 
            data.table()
        }
      }))
    
    # Bind back together, and add together
    mort_data <- rbind(mort_data, mort_data_ad1_area)[, .(events = sum(events)), by = c("area", "year", "sex", "age")][order(area, year, sex, age, events)]
  }
  
  # Define number of rows if data was square (i.e. a value for each uid, year, sex, age) 
  # Notice that not all uid's are necessarily included 
  square_rows <- 
    mort_data %>% 
    dplyr::select(-events) %>% 
    summarise_all(funs(length(unique(.)))) %>% 
    prod()
  
  pct_missing <- 
    dplyr::count(mort_data, area, year, sex, age) %>% 
    {100*(1 - (nrow(.) / square_rows))} %>% 
    round(1)
  
  message(paste0(pct_missing, "% of area-year-sex-age observations are missing from mortality data,\nwe assume this is due to no deaths in this area-year-sex-age observation"))
  
  
  # Find all stable areas, this is done by loading shapefile becuase some areas are missing from mortality if they never record any deaths
  areas <- 
    st_read(vr_pull_shp(paste0(country, "_", admin), "stable", "full")$shapefile, quiet = T) %>% 
    pull(uid) %>% 
    unique() %>% as.character() %>% as.numeric() %>% sort()
  
  # Fill missing area-year-sex-age observations with zero. 
  events <-
    expand.grid(area =  areas,
                year = unique(mort_data$year),
                sex  = unique(mort_data$sex),
                age  = unique(mort_data$age)) %>%
    data.table() %>%
    left_join(mort_data, by = c("area", "year", "sex", "age")) %>%
    replace_na(list(events = 0)) %>%
    arrange(area, year, sex, age) %>%
    data.table()
  
  # Make sure data is square
  assert_that(events %>% dplyr::select(-events) %>%
                summarise_all(funs(length(unique(.)))) %>%
                prod() == nrow(events))
  return(events)
}

vr_plot_raking_factors <- function(pop_rf, save_file){

  # Make sure graph is centered at 1
  limit <- max(c(abs(1-min(pop_rf$rf)), abs(1-max(pop_rf$rf)))) + 0.05
  limit <- c(1 - limit, 1 + limit)
  # Create plot 
  gg <- 
    ggplot(pop_rf) + 
    geom_density(aes(x = rf, color = factor(year)), adjust = 1.5) +
    theme_classic() + 
    lims(x = limit) + 
    labs(x = "Raking factor", y = "Density") + 
    scale_color_discrete("Year") + 
    facet_wrap(~ age, scales = "free_y")
  ggsave(filename = save_file, gg)
}

# Plot aggregated covariates 
vr_plot_aggregated_covs <- function(agg_covar, shapefile_path, 
                                    shapefile_field, save_file) {
  plot_covs <- function(data, cov) {
    colors <- c('#ffffe0','#ffe4ac','#ffc879','#ffa84c','#ff8725',
                '#ff5c03','#f12861','#cc117d','#a60383','#800080')
    color_values <- rescale(unique(c(seq(0, 5, length.out = 4), 
                                     seq(5, 10, length.out = 4), 
                                     seq(10, 25, length.out = 4))))
    ggplot(data) +
      geom_sf(aes(fill = get(cov))) +
      scale_fill_gradientn(colors = rev(colors),
                           values = color_values,
                           name = eval(cov)) +
      theme_void() +
      coord_sf(datum = NA) +
      facet_wrap(~ year)
  }
  
  polygon <- 
    st_read(shapefile_path, quiet = T) %>% 
    rename(area = shapefile_field) %>% 
    dplyr::mutate(area = as.numeric(as.character(area)))

  data <- 
    covar %>% 
    left_join(polygon, by = "area")
  
  covs <- names(covar) %>% setdiff(c("area", "year"))
  message("Saving")
  # Save to pdf 
  pdf(save_file, height = 10, width = 10)
  map(covs, ~ plot_covs(data, .x)) %>% 
    walk(print)
  dev.off()
  message("Done")
}

# Plot HIV mortality data 
vr_plot_mortality_data <- function(mort_data, save_file) {

  mor1 <-
    mort_data %>%
    group_by(year, age, sex) %>%
    dplyr::summarize(events = sum(events)) %>%
    ungroup() %>%
    mutate(sex_label = factor(sex, 1:2, c("Males", "Females"))) %>%
    ggplot(aes(x = year, y = events, color = sex_label)) +
    geom_point() +
    stat_smooth(geom = 'line', alpha = 0.5, se = FALSE) +
    facet_wrap( ~ age) +
    theme_bw() +
    labs(y = "HIV deaths by age")
  
  mor2 <-
    mort_data %>%
    group_by(year, age, sex) %>%
    dplyr::summarize(events = sum(events)) %>%
    ungroup() %>%
    mutate(sex_label = factor(sex, 1:2, c("Males", "Females"))) %>%
    ggplot(aes(x = age, y = events, color = sex_label)) +
    geom_point() +
    stat_smooth(geom = 'line', alpha = 0.5, se = FALSE) +
    facet_wrap( ~ year) +
    theme_bw() +
    labs(y = "HIV deaths by age")
  
  mor3 <-
    mort_data %>%
    group_by(year, age, sex) %>%
    dplyr::summarize(events = sum(events)) %>%
    ungroup() %>%
    mutate(sex_label = factor(sex, 1:2, c("Males", "Females"))) %>%
    ggplot(aes(x = year, y = events, color = age)) +
    geom_point() +
    stat_smooth(geom = 'line', alpha = 0.5, se = FALSE) +
    facet_wrap( ~ sex) +
    theme_bw() +
    labs(y = "HIV deaths by age")
  
  pdf(save_file, height = 10, width = 10)
  print(mor1)
  print(mor2)
  print(mor3)
  dev.off()
  message("Done plotting!")
  
}

# Function to prep u5m mortality data to be used in HIV small area estimation models for completeness assumption

#' @title vr_prep_hiv_mortality
#' @description: This function preps VR mortality data to be used in HIV small area estimation models. 
#' 
#'
#' @param country iso3 code for the country that you wish to prep HIV mortality data on
#' @param admin Admin level (for example "adm2")
#' @param years Years of population estimates 
#' @param shapefile_path Path to shapefile used 
#' @param agg_level This is the level of aggregation of prepped all cause mortality, either admin1 or area
#' 
#' @return Returns a data table that has the each admin unit by year and the number of u5m deaths after redistribution

vr_prep_all_cause_mortality <- function(country, 
                                        admin,
                                        years,
                                        agg_level = "area",
                                        cause = "all_cause") {
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(assertthat)
  summarize <- dplyr::summarize
  
  # Link between GBD age group and SAE ages
  age_link <- 
    get_age_metadata(age_group_set_id = 12, gbd_round_id = 6) %>% 
    dplyr::select(age_group_id, age_group_years_start) %>% 
    dplyr::rename(age = age_group_years_start) %>% 
    data.table()
  
  # Link between admin code and GBD location id
  core_repo <- paste0("<<<< FILEPATH REDACTED >>>>/lbd_core/")
  source(paste0(core_repo, "data_central/vr_functions.R"))
  
  # Link between admin2 code and GBD location ID 
  loc_id_to_area <-
    vr_pull_loc_meta(paste0(country, "_", admin)) %>%
    filter(level == as.numeric(str_extract(admin, "[0-9]+"))) %>%
    dplyr::select(area = uid, location_id) %>%
    arrange(area, location_id) %>% 
    unique()
  
  # Load country specific raw mortlity data (this may take some time)
  mort_raw <- vr_pull_cod(paste0(country,"_", admin), "redistribution")
 
  # Subset to correct years
  mort_raw <- mort_raw[year_id %in% years,]
  
  # See what percent of HIV deaths are mapped to country code or admin1 instead of admin2
  country_code <- get_gbd_locs(country, shapefile_version = "current", rake_subnational = F)[, location_id]
  
  pct_country <- 
    mort_raw %>%
    mutate(level = ifelse(location_id %in% country_code, 1, 2)) %>%
    group_by(level) %>%
    dplyr::summarize(deaths = sum(deaths)) %>%
    ungroup() %>%
    mutate(pct = round(100*deaths / sum(deaths), 3)) %>%
    data.table()
  
  admin1_code <- 
    vr_pull_loc_meta(paste0(country, "_", admin)) %>% 
    filter(level == 1) %>% 
    pull(location_id)
  
  pct_admin1 <-
    mort_raw %>%
    mutate(level = ifelse(location_id %in% admin1_code, 1, 2)) %>%
    group_by(level) %>%
    summarize(deaths = sum(deaths)) %>%
    ungroup() %>%
    mutate(pct = round(100*deaths / sum(deaths), 3)) %>%
    data.table()
  
 # Check how many deaths are assigned to admin1 or country 
  if (nrow(filter(mort_raw, location_id %in% country_code)) != 0) {
    message(paste0(pct_country[level == 1, pct], "% of deaths were mapped to country instead of admin 2. Currently those are being dropped"))
  }
  
  if (nrow(filter(mort_raw, location_id %in% admin1_code)) != 0) {
    message(paste0(pct_admin1[level == 1, pct], "% of deaths were mapped to admin 1 instead of admin 2. Currently those are being dropped"))
  }
  
  # Make sure all location ids in redistribution are included in location metadata
  # Sometimes matched to area 
  if (!are_equal(mort_raw$location_id %>% unique() %>% sort(), loc_id_to_area$location_id %>% unique() %>% sort())) {
    message(paste("Warning there differences between the location id's in the redistributed mortality data and the location metadata, look into this"))
    setdiff(mort_raw$location_id %>% unique(), loc_id_to_area$location_id %>% unique())
    setdiff(loc_id_to_area$location_id %>% unique(), mort_raw$location_id %>% unique())
  }
  
  
  # Get correct age grouping and area ID for HIV  ----------------------------------------
  # For now omit things mapped to admin 1
  # Take raw mortality data and collapse by year and location id
  if (agg_level == "area"){
    
    mort_data <- 
      mort_raw %>%
      filter(!location_id %in% admin1_code) %>% 
      left_join(age_link, by = "age_group_id") %>% # get SAE age groupings, five years up until 80+
      left_join(loc_id_to_area, by = "location_id") %>% # change to area
      dplyr::select(area, age, year = year_id, sex = sex_id, deaths) %>% 
      group_by(area, age, year, sex) %>% # collapse by uid, year, sex, age
      dplyr::summarize(deaths = sum(deaths)) %>% 
      ungroup() %>% 
      arrange(area, age, year, sex, deaths) %>% 
      data.table()
    
    # Add on admin-1 if required 
    if (nrow(mort_raw[location_id %in% admin1_code]) != 0) { 
      message("Assigning extra deaths to admin-1 to admin-2 proportional to deaths at that level")
      admin1_area_map <- 
        vr_pull_loc_meta(paste0(country, "_", admin)) %>% 
        filter(level == 2, year_end == max(year_end)) %>% # Gets most recent map to 
        dplyr::select(area = uid, admin2 = location_id, admin1 = location_parent_id) %>% 
        data.table()
      
      # Grab all deaths from raw data that are coded to admin-1
      mort_data_ad1 <-
        mort_raw %>% 
        filter(location_id %in% admin1_code) %>% 
        left_join(age_link, by = c("age_group_id")) %>% 
        mutate(age = case_when(age < 5 ~ 0,
                               between(age, 5, 75) ~ age,
                               age >= 80 ~ 80)) %>% 
        dplyr::select(admin1 = location_id, year = year_id, sex = sex_id, age, deaths) %>% 
        group_by(admin1, year, sex, age) %>% 
        dplyr::summarize(deaths = sum(deaths)) %>% 
        filter(deaths != 0) %>% # Filter out deaths that are zero
        data.table()
      
      # Redistribute deaths for each area-year-age-sex proportional to how many deaths are already coded there
      mort_data_ad1_area <-
        rbindlist(mclapply(split(mort_data_ad1, 1:nrow(mort_data_ad1)), mc.cores = 25, function(x) {
          area_specific <-
            mort_data[area %in% admin1_area_map[admin1 == x[["admin1"]], area] & year == x[["year"]] & sex == x[["sex"]] & age == x[["age"]]]
          if (nrow(area_specific[deaths != 0]) != 0) {
            area_specific <-
              area_specific %>%
              mutate(pct = deaths / sum(deaths)) %>%
              mutate(redistribute_deaths = x[["deaths"]]) %>%
              mutate(deaths = redistribute_deaths * pct) %>%
              dplyr::select(area, year, sex, age, deaths) %>% 
              data.table()
          } else {
            # If there were no deaths recorded in area-year-age-sex then split equally
            area_specific <-
              area_specific %>%
              mutate(deaths = x[["deaths"]] / nrow(area_specific)) %>% 
              dplyr::select(area, year, sex, age, deaths) %>% 
              data.table()
          }
        }))
      
      # Bind back together, and add together
      mort_data <- rbind(mort_data, mort_data_ad1_area)[, .(deaths = sum(deaths)), by = c("area", "year", "sex", "age")][order(area, year, sex, age, deaths)]
    }
    
    
  } else if (agg_level == "admin1") {
    message("Aggregation level is admin1")
    # Grab deaths mapped to admin1
    mort_admin1 <- merge(mort_raw[location_id %in% admin1_code],
                         age_link, by = "age_group_id")[, .(area = location_id, age, year = year_id, sex = sex_id, deaths)]
    
    # Assign deaths mapped to area to their respective admin1
    area_admin1_map <- 
      vr_pull_loc_meta(paste0(country, "_", admin)) %>% 
      filter(level == 2, year_end == max(year_end)) %>% # Gets most recent map to 
      dplyr::select(area = uid, admin1 = location_parent_id) %>% 
      data.table()
    
    mort_area_admin1 <- 
      mort_raw[!location_id %in% admin1_code, ] %>% 
      left_join(age_link, by = "age_group_id") %>% # get SAE age groupings, five years up until 80+
      left_join(loc_id_to_area, by = "location_id") %>% # change to area
      dplyr::select(area, age, year = year_id, sex = sex_id, deaths) %>% 
      group_by(area, age, year, sex) %>% # collapse by uid, year, sex, age
      dplyr::summarize(deaths = sum(deaths)) %>% 
      ungroup() %>% 
      arrange(area, age, year, sex, deaths) %>% 
      data.table()
    
    mort_area_admin1 <- 
      mort_area_admin1 %>% 
      filter(!is.na(area)) %>% 
      left_join(area_admin1_map, by = "area") %>% 
      dplyr::select(area = admin1, age, year, sex, deaths) %>% 
      data.table()
      
    mort_data <- rbind(mort_admin1, mort_area_admin1)[, .(deaths = sum(deaths)), .(area, age, year, sex)]
    assert_that(mort_data %>% dplyr::select(-deaths) %>%
      summarise_all(funs(length(unique(.)))) %>%
      prod() == nrow(mort_data))
    return(mort_data)
  } else {
    stop("Agg level only supports area and admin1 at this time")
  }
  
  # This is because of some areas not being in location metadata table (looking into)
  message(round(100*sum(mort_data[is.na(area), deaths]) / sum(mort_data[, deaths]), 3), "% of deaths being dropped becuase of missing area")
  mort_data <- mort_data[!is.na(area)]
  
  # Define number of rows if data was square (i.e. a value for each uid, year, sex, age) 
  # Notice that not all uid's are necessarily included 
  square_rows <- 
    mort_data %>% 
    dplyr::select(-deaths) %>% 
    summarise_all(funs(length(unique(.)))) %>% 
    prod()
  
  pct_missing <- 
    dplyr::count(mort_data, area, age, year, sex) %>% 
    {100*(1 - (nrow(.) / square_rows))} %>% 
    round(1)
  
  message(paste0(pct_missing, "% of area-year-sex-age observations are missing from mortality data, we assume this is due to no deaths in this area-year-sex-age observation"))
  
  # Find all stable areas, this is done by loading shapefile becuase some areas are missing from mortality if they never record any deaths
  areas <- 
    st_read(vr_pull_shp(paste0(country, "_", admin), "stable", "full")$shapefile, quiet = T) %>% 
    pull(uid) %>% 
    unique() %>% as.character() %>% as.numeric() %>% sort()
  
  # Fill missing area-year-sex-age observations with zero. 
  vr_deaths <-
    expand.grid(area =  areas,
                age  = unique(mort_data$age), 
                year = unique(mort_data$year),
                sex  = unique(mort_data$sex)) %>% 
    left_join(mort_data, by = c("area", "age", "year", "sex")) %>%
    replace_na(list(deaths = 0)) %>%
    arrange(area, age, year, sex) %>%
    data.table()
  
  # Make sure data is square
  assert_that(vr_deaths %>% dplyr::select(-deaths) %>%
                summarise_all(funs(length(unique(.)))) %>%
                prod() == nrow(vr_deaths))
  return(vr_deaths)
}

vr_plot_raking_factors <- function(pop_rf, save_file){
  
  # Make sure graph is centered at 1
  limit <- max(c(abs(1-min(pop_rf$rf)), abs(1-max(pop_rf$rf)))) + 0.05
  limit <- c(1 - limit, 1 + limit)
  # Create plot 
  gg <- 
    ggplot(pop_rf) + 
    geom_density(aes(x = rf, color = factor(year)), adjust = 1.5) +
    theme_classic() + 
    lims(x = limit) + 
    labs(x = "Raking factor", y = "Density") + 
    scale_color_discrete("Year") + 
    facet_wrap(~ age, scales = "free_y")
  ggsave(filename = save_file, gg)
}

# Plot aggregated covariates 
vr_plot_aggregated_covs <- function(agg_covar, shapefile_path, 
                                    shapefile_field, save_file) {
  plot_covs <- function(data, cov) {
    colors <- c('#ffffe0','#ffe4ac','#ffc879','#ffa84c','#ff8725',
                '#ff5c03','#f12861','#cc117d','#a60383','#800080')
    color_values <- rescale(unique(c(seq(0, 5, length.out = 4), 
                                     seq(5, 10, length.out = 4), 
                                     seq(10, 25, length.out = 4))))
    ggplot(data) +
      geom_sf(aes(fill = get(cov)), color = NA) +
      scale_fill_gradientn(colors = rev(colors),
                           values = color_values,
                           name = eval(cov)) +
      theme_void() +
      coord_sf(datum = NA) +
      facet_wrap(~ year)
  }
  
  polygon <- 
    st_read(shapefile_path, quiet = T) %>% 
    dplyr::rename(area = shapefile_field) %>% 
    dplyr::mutate(area = as.numeric(as.character(area)))
  
  data <- 
    polygon %>% 
    left_join(covar)
  
  covs <- names(covar) %>% setdiff(c("area", "year"))
  message("Saving")
  # Save to pdf 
  pdf(save_file, height = 10, width = 10)
  map(covs, ~ plot_covs(data, .x)) %>% 
    walk(print)
  dev.off()
  message("Done")
}
