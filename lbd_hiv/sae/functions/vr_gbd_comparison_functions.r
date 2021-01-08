## Function that plots comparison between model fits unraked and GBD mortality
# This plot sees how the model fits the data and compares to GBD


admin_0_gbd_compare <- function(country, 
                                run_date,
                                rake = T,
                                core_repo = paste0("<<<< FILEPATH REDACTED >>>>/lbd_core/"), 
                                gbd = 2017) {
  
  
  library(data.table)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  library(tidyr)
  library(sf)

  select    <- dplyr::select
  rename    <- dplyr::rename
  summarize <- dplyr::summarize
  source(paste0(core_repo, '/mbg_central/post_estimation_functions.R'))
  source(paste0(core_repo, '/mbg_central/shapefile_functions.R'))
  source(paste0(core_repo, "mbg_central/misc_functions.R"))
  
  # Source VR functions 
  function_files = list.files(paste0("<<<< FILEPATH REDACTED >>>>/lbd_hiv/sae/functions/"), full.names = T)
  sapply(function_files, source)
  
  # Load years and ages 
  model_dir <- paste0("<<<< FILEPATH REDACTED >>>>")
  source(paste0(core_repo, "/sae_central/settings.r"))
  get_settings(model_dir)
  
  ## Load final model estimates unraked ---------------------------------------------------
  load(paste0(model_dir,"/est_all.rdata"))
  est <- est[sex < 3 & level == "national",]
  est[, sex_label := factor(sex, 2:1, c("Females", "Males"))]
  setnames(est, "mean", "sae")

  # Grab GBD location id(s)
  loc_ids <- get_gbd_locs(country, shapefile_version = "current", rake_subnational = T)$location_id
  gbd_input <- load_gbd_estimates(loc_ids, years, cause_id = "298", include_data = T, gbd = gbd)
  
  gbd_data <- gbd_input[["data"]]
  gbd_estimates <- gbd_input[["estimates"]]
  
  # Load model data ----------------------------------------------------------------
  age_std <- fread(age_std_file)
  
  load(paste0(model_dir, '/temp_dir/data.rdata'))
  data[, year := year + min(years)]
  data[, age := ages[age + 1]]
  data_unstandard <- data[, list(sae = sum(events) / sum(pop)), by = 'year,sex']
  data_unstandard[, age := 98]
  data <- data[, list(sae = sum(events) / sum(pop)), by = 'year,sex,age']
  
  # Get age standard of data 
  data_standard <- 
    data %>% 
    left_join(age_std) %>% 
    group_by(year, sex) %>% 
    dplyr::summarize(sae = weighted.mean(sae, wt)) %>% 
    ungroup() %>% 
    mutate(age = 99) %>% 
    data.table()
  
  # Join data together
  data <- rbind(data, data_unstandard, data_standard)
  
  data <-
    data %>%
    left_join(gbd_data, by = c("year", "age", "sex")) %>%
    mutate_at(vars(gbd, sae), funs(.*100000)) %>%
    gather(key = "Source", value = "data", sae, gbd) %>%
    mutate(Source = factor(Source, c("sae", "gbd"), c("SAE", "GBD"))) %>%
    mutate(sex_label = factor(sex, 1:2, c("Males", "Females")))
     
  # Make plot of GBD trend line compared to LBD
  gg1 <- 
    gbd_estimates %>% 
    left_join(est, c("age", "year", "sex")) %>% 
    mutate(sex_label = factor(sex, 1:2, c("Males", "Females"))) %>%
    mutate_at(vars(gbd, sae), funs(.*100000)) %>% 
    gather(key = "Source", value = "mortality", sae, gbd) %>%
    mutate(Source = factor(Source, c("sae", "gbd"), c("SAE", "GBD"))) %>% 
    ggplot() + 
    geom_line(aes(x = year, y = mortality, color = sex_label, linetype = Source)) + 
    geom_point(data = data, aes(x = year, y = data, color = sex_label, shape = Source)) +
    scale_color_discrete(name = "Sex") +
    labs(y = "HIV deaths per 100,000 people", x = "Year") +
    theme_light() + 
    facet_wrap(~ age, scales = "free")
  
  # Make unraked estimates with data at nationa level, with uncertainty included
  gg_unraked <- 
    est %>% 
    mutate_at(vars(sae, lb, ub), funs(.*100000)) %>% 
    ggplot() + 
    geom_line(aes(x = year, y = sae, color = sex_label)) + 
    geom_ribbon(aes(x = year, ymin = lb, ymax = ub, fill = sex_label), alpha = 0.2) + 
    geom_point(data = filter(data, Source == "SAE"), aes(x = year, y = data, color = sex_label), shape = 1) +
    scale_color_discrete(name = "Sex") +
    scale_fill_discrete(name = "Sex") +
    labs(y = "HIV deaths per 100,000 people", x = "Year") +
    theme_light() + 
    facet_wrap(~ age, scales = "free")
  
  if (rake) {
    # Load raked estimates
    load(paste0(model_dir,"/est_all_raked.rdata"))
    est <- est[sex < 3 & level == "national",]
    est[, sex_label := factor(sex, 2:1, c("Females", "Males"))]
    setnames(est, "mean", "sae")

    # Save raked estimates
    gg_raked <- 
      est %>% 
      mutate(sex_label = factor(sex, 1:2, c("Males", "Females"))) %>% 
      mutate_at(vars(sae, lb, ub), funs(.*100000)) %>% 
      ggplot(aes(x = year, y = sae)) + 
      geom_point(data = filter(data, Source == "SAE"), aes(x = year, y = data, color = sex_label), shape = 1) +
      geom_line(aes(color = sex_label)) + 
      geom_ribbon(aes(ymin = lb, ymax = ub, fill = sex_label), alpha = 0.2) + 
      scale_color_discrete(name = "Sex") +
      scale_fill_discrete(name = "Sex") +
      labs(y = "HIV deaths per 100,000 people", x = "Year") +
      theme_light() + 
      facet_wrap(~ age, scales = "free")
    
  }
  
  # Save graphs 
  dir.create(paste0(model_dir, "/rake_plots/"), showWarnings = FALSE)
  pdf(file = paste0(model_dir, "/rake_plots/gbd_compare_data_", gbd, ".pdf"), height = 8, width = 12)
  print(gg1)
  print(gg_unraked)
  if (rake) print(gg_raked)
  dev.off()

  message(paste0("Plot saved to ", paste0(model_dir, "gbd_compare_data.pdf")))
}

admin_0_completeness_compare <- function(country,
                                         run_dates,
                                         run_date_labels = c("standard", "completeness"),
                                         gbd = 2019) {
  
  
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  library(tidyr)
  library(sf)
  
  core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>/lbd_core/")
  source('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/post_estimation_functions.R')
  source('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/shapefile_functions.R')
  source(paste0(core_repo, "mbg_central/misc_functions.R"))
  source(paste0(core_repo, '/mbg_central/setup.R'))
  package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv")) %>% setdiff(("plyr"))
  mbg_setup(package_list = package_list, repos = core_repo)
  
  select    <- dplyr::select
  rename    <- dplyr::rename
  summarize <- dplyr::summarize

  # Load years and ages 
  raw_est <- 
    rbindlist(lapply(1:length(run_dates), function(i) {
      # Load years and ages 
      model_dir <- paste0("<<<< FILEPATH REDACTED >>>>")
      source(paste0(core_repo, "/sae_central/settings.r"))
      get_settings(model_dir)
      
      ## Load final model estimates unraked ---------------------------------------------------
      load(paste0(model_dir,"/est_all.rdata"))
      est <- est[sex %in% c(1,2) & level == "national",]
      est[, run := run_date_labels[i]]
    }))
  est <- raw_est[, .(year, age, sex, mean, run)]
  
  # Grab GBD location id(s)
  loc_ids <- get_gbd_locs(country, shapefile_version = "current", rake_subnational = T)$location_id
  gbd_input <- load_gbd_estimates(loc_ids, years, cause_id = "298", include_data = T, gbd = gbd)
  
  gbd_data <- gbd_input[["data"]]
  gbd_estimates <- gbd_input[["estimates"]]
  
  gbd_estimates <- 
    gbd_estimates %>% 
    select(age, year, sex, mean = gbd) %>% 
    mutate(run = "gbd") %>% 
    data.table()
  
  # Bind GBD 
  est <- rbind(est, gbd_estimates)
  
  # Load model data ----------------------------------------------------------------
  model_dir <- paste0("<<<< FILEPATH REDACTED >>>>")
  source(paste0(core_repo, "/sae_central/settings.r"))
  get_settings(model_dir)
  
  
  age_std <- fread(age_std_file)
  load(paste0(model_dir, '/temp_dir/data.rdata'))
  data[, year := year + min(years)]
  data[, age := ages[age + 1]]
  data_unstandard_copy <- copy(data)
  data_unstandard <- data[, list(sae = sum(events) / sum(pop)), by = 'year,sex']
  data_unstandard[, age := 98]
  data <- data[, list(sae = sum(events) / sum(pop)), by = 'year,sex,age']
  
  # Get age standard of data 
  data_standard <- 
    data %>% 
    left_join(age_std) %>% 
    group_by(year, sex) %>% 
    dplyr::summarize(sae = weighted.mean(sae, wt)) %>% 
    ungroup() %>% 
    mutate(age = 99) %>% 
    data.table()
  
  # Join data together
  data <- rbind(data, data_unstandard, data_standard)
  
  # Bind together gbd data
  data <-
    data %>%
    left_join(gbd_data, by = c("year", "age", "sex")) %>%
    mutate_at(vars(gbd, sae), funs(.*100000)) %>%
    gather(key = "Source", value = "data", sae, gbd) %>%
    mutate(Source = factor(Source, c("sae", "gbd"), c("SAE", "GBD"))) %>%
    mutate(sex_label = factor(sex, 1:2, c("Males", "Females")))
  
 gg1 <- 
  est %>%
  mutate(sex_label = factor(sex, 1:2, c("Males", "Females"))) %>% 
  mutate_at(vars(mean), funs(.*100000)) %>% 
  mutate(run = factor(run, levels = c("completeness", "standard", "gbd"))) %>% 
  ggplot() + 
  geom_line(aes(x = year, y = mean, color = sex_label, linetype = run)) + 
  geom_point(data = data, aes(x = year, y = data, color = sex_label, shape = Source)) +
  scale_color_discrete(name = "Sex") +
  labs(y = "HIV deaths per 100,000 people", x = "Year") +
  theme_light() + 
  facet_wrap(~ age, scales = "free")
 
 ribbon_data <- 
   copy(raw_est[run == "completeness"]) %>% 
   mutate_at(vars(mean, lb, ub), funs(.*100000)) %>% 
   mutate(sex_label := factor(sex, 1:2, c("Males", "Females")))
 
 gg_unraked <- 
   est %>%
   mutate(sex_label = factor(sex, 1:2, c("Males", "Females"))) %>% 
   mutate_at(vars(mean), funs(.*100000)) %>% 
   filter(run %in% c("completeness", "gbd")) %>% 
   ggplot() + 
   geom_line(aes(x = year, y = mean, color = sex_label, linetype = run)) + 
   geom_point(data = data, aes(x = year, y = data, color = sex_label, shape = Source)) +
   geom_ribbon(data = ribbon_data,
               aes(x = year, ymin = lb, ymax = ub, fill = sex_label), alpha = 0.2) + 
   scale_color_discrete(name = "Sex") +
   labs(y = "HIV deaths per 100,000 people", x = "Year") +
   theme_light() + 
   facet_wrap(~ age, scales = "free")
  
  # Save graphs 
  dir.create(paste0(model_dir, "/completeness_plots/"), showWarnings = FALSE)
  pdf(file = paste0(model_dir, "/completeness_plots/gbd_completeness_comparison_", gbd, ".pdf"), height = 8, width = 12)
  print(gg1)
  print(gg_unraked)
  dev.off()
  
  message(paste0("Plot saved to ", paste0(model_dir, "/completeness_plots/gbd_completeness_comparison_", gbd, ".pdf")))
}


compare_raking_factors <- function(country,
                                   run_dates,
                                   run_date_labels = c("standard", "completeness"),
                                   return_plot = F,
                                   out_dir =  paste0("<<<< FILEPATH REDACTED >>>>")){
  
  int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0]
  library(tidybayes, lib.loc = "<<<< FILEPATH REDACTED >>>>")
  
  rf <- 
    rbindlist(lapply(1:length(run_dates), function(i) {
      # Load years and ages 
      model_dir <- paste0("<<<< FILEPATH REDACTED >>>>")
      source(paste0(core_repo, "/sae_central/settings.r"))
      get_settings(model_dir)
      
      rf <-
        rbindlist(lapply(years, function(year) {
          rbindlist(lapply(sexes, function(sex) {
            fread(paste0(model_dir, "/temp_dir/rf_", year, "_", sex, ".csv"))
          }))
        }))
     
      rf[, run := run_date_labels[i]]
    }))
  rf <- rf[, .(year, age, sex, rf, run)]
  
  rf <- 
    rf %>% 
    mutate(age = case_when(age == 0  ~ "<5",
                           age == 5  ~ "5-9",
                           age == 10 ~ "10-14", 
                           age == 15 ~ "15-19",
                           age == 20 ~ "20-24",
                           age == 25 ~ "25-29",
                           age == 30 ~ "30-34",
                           age == 35 ~ "35-39",
                           age == 40 ~ "40-44",
                           age == 45 ~ "45-49",
                           age == 50 ~ "50-54",
                           age == 55 ~ "55-59",
                           age == 60 ~ "60-64",
                           age == 65 ~ "65-69",
                           age == 70 ~ "70-74",
                           age == 75 ~ "75-79",
                           age == 80 ~ "80+")) %>% 
    mutate(age = factor(age, levels = c("<5", "5-9", "10-14", "15-19", "20-24",
                                        "25-29", "30-34", "35-39", "40-44", "45-49",
                                        "50-54", "55-59", "60-64", "65-69", "70-74", 
                                        "75-79", "80+")))
  
  gg_rf <- 
    ggplot(rf) + 
    geom_abline(intercept = 1, slope = 0, linetype = "dashed", size = 0.2, color = "black") + 
    stat_halfeye(aes(x = run, y = rf), side = "left", .width = c(0.95), size = 1) + 
    labs(x = NULL, y = "Raking factor") + 
    facet_wrap(~ age, scales = "free", ncol = 4) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black")) + 
    coord_cartesian(ylim = c(0, 3))
  
  if (return_plot) return(gg_rf)
  
  dir.create(out_dir, showWarnings = FALSE)
  pdf(file = paste0(out_dir, "/rf_comparison.pdf"), height = 11, width = 8.8)
  print(gg_rf)
  dev.off()
}

compare_gbd_runs <- function(country, out_dir) {
  
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  library(tidyr)
  library(sf)
  
  core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>/lbd_core/")
  source('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/post_estimation_functions.R')
  source('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/shapefile_functions.R')
  source(paste0(core_repo, "mbg_central/misc_functions.R"))
  source(paste0(core_repo, '/mbg_central/setup.R'))
  package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv")) %>% setdiff(("plyr"))
  mbg_setup(package_list = package_list, repos = core_repo)
  
  select    <- dplyr::select
  rename    <- dplyr::rename
  summarize <- dplyr::summarize
  
  # Grab GBD location id(s)
  loc_ids <- get_gbd_locs(country, shapefile_version = "current", rake_subnational = T)$location_id
  gbd_input_2017     <- load_gbd_estimates(loc_ids, years, cause_id = "298", include_data = T, gbd = 2017)
  gbd_data_2017      <- data.table(gbd_input_2017 [["data"]])
  gbd_estimates_2017 <- data.table(gbd_input_2017 [["estimates"]])
  
  gbd_input_2019     <- load_gbd_estimates(loc_ids, years, cause_id = "298", include_data = T, gbd = 2019)
  gbd_data_2019      <- data.table(gbd_input_2019[["data"]])
  gbd_estimates_2019 <- data.table(gbd_input_2019[["estimates"]])
  
  gbd_estimates_2017[, run := "GBD_2017"]
  gbd_estimates_2019[, run := "GBD_2019"]
  gbd_est <- 
    rbind(gbd_estimates_2017, gbd_estimates_2019) %>% 
    mutate(sex_label = factor(sex, 1:2, c("Males", "Females"))) %>% 
    mutate_at(vars(gbd), funs(.*100000)) %>% 
    mutate(run = factor(run, levels = c("GBD_2017", "GBD_2019")))
  
  gbd_data_2017[, run := "GBD_2017"]
  gbd_data_2019[, run := "GBD_2019"]
  gbd_data <- 
    rbind(gbd_data_2017, gbd_data_2019) %>% 
    mutate(run = factor(run, levels = c("GBD_2017", "GBD_2019"))) %>% 
    mutate(sex_label = factor(sex, 1:2, c("Males", "Females"))) %>% 
    mutate_at(vars(gbd), funs(.*100000))
  
  gg1 <- 
    ggplot(gbd_est) + 
    geom_line(aes(x = year, y = gbd, color = sex_label, linetype = run)) + 
    geom_point(data = gbd_data, aes(x = year, y = gbd, color = sex_label, shape = run)) +
    scale_color_discrete(name = "Sex") +
    labs(y = "HIV deaths per 100,000 people", x = "Year") +
    theme_light() + 
    facet_wrap(~ age, scales = "free")
  
  pdf(file = paste0(out_dir, "/GBD_comparison.pdf"), height = 8, width = 12)
  print(gg1)
  dev.off()
  
}

#### Used for loading GBD estimates and for finding K in completeness functions


# Loads estimates by age, year, sex of gbd mortality in correct age bins
# Default cause id is HIV
# if include data is true, append the input data to the gbd estimates

load_gbd_estimates <- function(loc_ids,
                               years,
                               metric_id = 3,
                               cause_id = "298",
                               include_data = F,
                               gbd = 2019) {
  
  library(dplyr)
  
  source("<<<< FILEPATH REDACTED >>>>/get_age_metadata.R")
  age_link <- 
    get_age_metadata(age_group_set_id = 12, gbd_round_id = 6) %>% 
    filter(between(age_group_years_start, 5, 75)) %>% 
    dplyr::select(age_group_id, age_group_years_start) %>% 
    dplyr::rename(age = age_group_years_start)
  
  # Grab GBD ids for under 5 and over 80 year olds
  under5 <- data.table(age_group_id = 1, age = 0)
  over80 <- data.table(age_group_id = 21, age = 80)
  age_link <- rbind(under5, age_link, over80)
  
  # Get outputs 
  source("<<<< FILEPATH REDACTED >>>>/get_outputs.R")
  source("<<<< FILEPATH REDACTED >>>>/get_ids.R")  
  #get_ids("age_group")
  
  # Not sure which output I have for GBD
  if (gbd == 2019) {
    gbd_output <- 
      get_outputs("cause",
                  year_id = years,
                  gbd_round_id = 6,
                  location_id = loc_ids,
                  sex_id = c(1, 2),
                  decomp_step = 'step5',
                  version = "latest", 
                  age_group_id = age_link$age_group_id,
                  cause_id = cause_id,
                  metric_id = metric_id) %>% 
      left_join(age_link, by = "age_group_id") %>% 
      dplyr::select(age, year = year_id, sex = sex_id, loc_id = location_id, gbd = val)
    
    age_link <- data.table(age = c(98, 99), age_group_id = c(22, 27))
    gbd_output_all_age <- 
      get_outputs("cause",
                  year_id = years,
                  gbd_round_id = 6,
                  location_id = loc_ids,
                  sex_id = c(1, 2),
                  decomp_step = 'step5',
                  version = "latest", 
                  age_group_id = age_link$age_group_id,
                  cause_id = cause_id,
                  metric_id = metric_id) %>% 
      left_join(age_link, by = "age_group_id") %>% 
      dplyr::select(age, year = year_id, sex = sex_id, loc_id = location_id, gbd = val)
    
    gbd_output <- rbind(gbd_output, gbd_output_all_age)
  } else if (gbd == 2017) {
    gbd_output <- 
      get_outputs("cause",
                  year_id = years,
                  gbd_round_id = 5,
                  location_id = loc_ids,
                  sex_id = c(1, 2),
                  version = "latest", 
                  age_group_id = age_link$age_group_id,
                  cause_id = cause_id,
                  metric_id = metric_id) %>% 
      left_join(age_link, by = "age_group_id") %>% 
      dplyr::select(age, year = year_id, sex = sex_id, loc_id = location_id, gbd = val)
    
    age_link <- data.table(age = c(98, 99), age_group_id = c(22, 27))
    gbd_output_all_age <- 
      get_outputs("cause",
                  year_id = years,
                  gbd_round_id = 5,
                  location_id = loc_ids,
                  sex_id = c(1, 2),
                  version = "latest", 
                  age_group_id = age_link$age_group_id,
                  cause_id = cause_id,
                  metric_id = metric_id) %>% 
      left_join(age_link, by = "age_group_id") %>% 
      dplyr::select(age, year = year_id, sex = sex_id, loc_id = location_id, gbd = val)
    
    gbd_output <- rbind(gbd_output, gbd_output_all_age)
  } else {
    stop("GBD must be either 2017 or 2019 to rake to")
  }
  
  # Return output and data if requested
  if (!include_data) {
    return(gbd_output)
  } else {
    
    # Grab HIV input data, need different age mapping for get_cod_data since it has new age groups
    age_link <-
      get_age_metadata(age_group_set_id = 12, gbd_round_id = 6) %>%
      dplyr::select(age_group_id, age_group_years_start) %>%
      dplyr::rename(age = age_group_years_start)
    
    # Problems with gbd input 
    source("<<<< FILEPATH REDACTED >>>>/get_cod_data.R")
    if (gbd == 2019)
      gbd_input  <-
        get_cod_data(
          cause_id = cause_id,
          cause_set_id = 2,
          gbd_round_id = 6,
          decomp_step = 'step4',
          year_id = years,
          location_id = loc_ids,
          sex_id = c(1, 2)
        )
    else {
      gbd_input  <-
        get_cod_data(
          cause_id = cause_id,
          cause_set_id = 2,
          gbd_round_id = 5,
          year_id = years,
          location_id = loc_ids,
          sex_id = c(1, 2)
        )
    }
    
    gbd_input <- 
      gbd_input %>% 
      left_join(age_link, by = "age_group_id") %>% # only keep most detailed age groups
      filter(!is.na(age)) %>% 
      dplyr::select(age, year, sex, loc_id = location_id, gbd = rate, pop) %>% 
      mutate(age_group = case_when(age < 5  ~  0,
                                   age > 75 ~ 80,
                                   TRUE     ~ age)) %>% 
      group_by(age_group, year, sex, loc_id) %>% 
      summarize(gbd = weighted.mean(gbd, pop)) %>% 
      ungroup() %>% 
      dplyr::select(age = age_group, year, sex, loc_id, gbd)
    
    # Add all age and age standardized
    age_link <- data.table(age = c(98, 99), age_group_id = c(22, 27))
    if (gbd == 2019) {
      gbd_all_ages <-
        get_cod_data(
          cause_id = cause_id,
          cause_set_id = 2,
          gbd_round_id = 6,
          decomp_step = 'step4',
          year_id = years,
          location_id = loc_ids,
          sex_id = c(1, 2)
        )
    } else {
      gbd_all_ages <-
        get_cod_data(
          cause_id = cause_id,
          cause_set_id = 2,
          gbd_round_id = 5,
          year_id = years,
          location_id = loc_ids,
          sex_id = c(1, 2)
        )
    }
    
    gbd_all_ages <- 
      gbd_all_ages %>% 
      left_join(age_link, by = "age_group_id") %>% # only keep most detailed age groups
      filter(!is.na(age), !is.na(rate)) %>% 
      dplyr::select(age, year, sex, loc_id = location_id, gbd = rate)
    
    age_std <- fread("<<<< FILEPATH REDACTED >>>>")
    gbd_input_age_std <- 
      gbd_input %>% 
      left_join(age_std) %>% 
      group_by(year, sex, loc_id) %>% 
      dplyr::summarize(gbd = weighted.mean(gbd, wt)) %>% 
      ungroup() %>% 
      mutate(age = 99) %>% 
      data.table()
    
    gbd_input <- rbind(gbd_input, gbd_all_ages, gbd_input_age_std)
    
    return(list(data = gbd_input, estimates = gbd_output))
  }
}


