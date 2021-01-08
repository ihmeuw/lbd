### Make time trend plots with data ------------------------------------------------------------------------------


hiv_sae_time_trend <- function(country, 
                               admin = "adm2",
                               shapefile_field = "uid",
                               run_date,
                               core_repo = paste0("<<<< FILEPATH REDACTED >>>>/lbd_core/"), 
                               cores = 1) {
  
  # Load libraries
  library(grid)
  library(gridExtra)
  library(parallel)
  library(scales)
  library(ggrepel)
  library(stringr)
  library(sf)
  library(dplyr)
  library(nngeo, lib.loc = "<<<< FILEPATH REDACTED >>>>")
  
  # source functions
  source(paste0(core_repo, "data_central/vr_functions.R"))
  source(paste0(core_repo, "data_central/sae_functions.R"))
  
  # Functions
  select <- dplyr::select
  
  country_name <- 
    st_read(get_admin_shapefile(admin_level = 0, raking = F, version = "current")) %>% 
    filter(ADM0_CODE == get_adm0_codes(country)) %>% 
    pull(ADM0_NAME)
  
  # Define directories
  model_dir <- paste0("<<<< FILEPATH REDACTED >>>>")
  temp_dir <- paste0(model_dir, "/temp_dir/")

  # Check settings 
  source(paste0(core_repo, "/sae_central/settings.r"))
  get_settings(model_dir)

  # Load and simplify shapefiles ---------------------------------------------
  # Load admin0 shapefile
  ad0_shape_simple <- 
    st_read(get_admin_shapefile(admin_level = 0, raking = F, version = "current")) %>% 
    filter(ADM0_CODE == get_adm0_codes(country)) %>% 
    st_simplify(preserveTopology = T) %>% 
    dplyr::select(ADM0_CODE)
  
  # Get admin1 names
  admin1_names <- 
    vr_pull_loc_meta(paste0(country, "_adm2")) %>% 
    filter(level == 1) %>% 
    dplyr::select(admin1 = location_id, admin1_name = adm_name) %>% 
    data.table()
    
  # Link between area and admin1 code
  load(geoagg_files[["admin1"]])
  area_admin1_link <- unique(weights[, .(area, admin1)])

  ad1_shape_simple <- 
    vr_pull_shp(source = paste0(country, "_adm2"), 
                type = "annual", 
                year = "last", 
                admin_level = 1,
                return_object = "sf",
                link = F)$shapefile %>% 
    transmute(admin1 = GAUL_CODE) %>% 
    left_join(admin1_names, by = "admin1") %>% 
    nngeo::st_remove_holes() %>% 
    st_simplify(preserveTopology = T) 
 
  # Shapefile for areas
  area_shape_simple <- 
    st_read(vr_pull_shp(paste0(country, "_", admin), "stable", "full")$shapefile, quiet = T) %>% 
    select(area = uid) %>% 
    arrange(area) %>% 
    st_simplify(preserveTopology = T)
  
  # Grab admin population for calculating percentile population for reference
  load(pop_file)
  pop <- get(load(pop_file))
  
  load(geoagg_files[["admin1"]])
  
  pop_admin1_pct <- 
    copy(weights) %>% 
    group_by(admin1) %>% 
    summarize(pop = sum(pop)) %>% 
    ungroup() %>% 
    mutate(percentile = 100*rank(pop) / nrow(.)) %>% 
    arrange(percentile) %>% data.table()
  
  # Calculate mean population for areas
  pop_area_pct <-
    pop %>%
    group_by(area) %>%
    dplyr::summarize(pop = sum(pop)) %>%
    ungroup() %>%
    mutate(percentile = 100*rank(pop) / nrow(.)) %>%
    arrange(percentile) %>% data.table()
  
  
  # Create directory to save plots
  dir.create(paste0(model_dir, "/time_trend/"))
  
  # Plot all admin 1 time trends -----------------------------------------------------------------------
  
  # Load age standardization for age standardized data inputs
  age_wts <- fread(age_std_file)
  
  # Load input data 
  load(paste0(temp_dir, '/data.rdata'))
  data[, year := year + min(years)]
  data[, age := ages[age + 1]]
  data[, events := 100000*events / pop]
  
  # Get data for each admin1 weighted by population and age
  admin1_data <- 
    data %>% 
    left_join(area_admin1_link, by = "area") %>% 
    group_by(admin1, year, sex, age) %>% 
    summarize(events = weighted.mean(events, pop), # pop-weighted death rate per admin1 for each age group
              pop = sum(pop)) %>% 
    ungroup()
  
  # Get unstandard data (population weighted by ages)
  admin1_data_unstandard <- 
    admin1_data %>% 
    group_by(admin1, year, sex) %>% 
    summarize(events = weighted.mean(events, pop), # age-standardized death rate
              pop = sum(pop)) %>% 
    ungroup() %>% 
    mutate(age = 98)
  
  # Grab admin1 age standardized
  admin1_data_standard <- 
    admin1_data %>% 
    left_join(age_wts, by = "age") %>% 
    group_by(admin1, year, sex) %>% 
    summarize(events = weighted.mean(events, wt), # age-standardized death rate
              pop = sum(pop)) %>% 
    ungroup() %>% 
    mutate(age = 99)
  
  # Merget together to get age standard and age standardized mortality rates 
  admin1_data_complete <- 
    rbind(admin1_data, admin1_data_unstandard, admin1_data_standard) %>% 
    arrange(admin1, year, sex, age) %>% 
    data.table()
  
  # Load data estimates ---------------------------------------
  load(paste0(model_dir, "/est_all.rdata"))
  est <- est[sex < 3,]
  #est[, sex_label := factor(sex, 2:1, c("Females", "Males"))]
  
  # Estimates for admin1
  admin1_est <- 
    est[level == "admin1"] %>% 
    rename(admin1 = area) %>% 
    left_join(admin1_names, by = "admin1") %>% 
    mutate_at(vars(mean, lb, ub), funs(.*100000)) %>% 
    data.table()
  
  # Join together 
  admin1_est_data <- merge(admin1_est, admin1_data_complete, by = c("admin1", "year", "age", "sex"))
    
  # Make all admin1 graph  -------------------------------------------------------
  message("Saving all admin1 graph")
  int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0] 
  ad1_time_trend <-
    admin1_est_data %>% 
    mutate(sex_label = factor(sex, 1:2, c("Males", "Females"))) %>% 
    filter(age == 99) %>% 
    ggplot() + 
    geom_line(aes(x = year, y = mean, color = sex_label)) + 
    geom_ribbon(aes(x = year, ymin = lb, ymax = ub, fill = sex_label), alpha = 0.1, color = NA) + 
    geom_point(aes(x = year, y = events, fill = sex_label, size = pop), shape = 21, alpha = 0.5) +
    facet_wrap(~ admin1_name) +
    labs(title = paste("Age-standardized HIV deaths in", country_name, "by Admin-1"),
         x = "Year",
         y = "HIV deaths per 100,000 people") +
    scale_color_discrete(name = "Sex") +
    scale_fill_discrete(name = "Sex") + 
    scale_x_continuous(breaks = int_breaks) +
    scale_size_continuous(name = "Population", range = c(0.5,3)) + 
    theme_bw(base_size = 16) + 
    theme(strip.background = element_blank(),
          plot.caption = element_text(hjust = 0.5),
          plot.title = element_text(size = 15, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 7),
          legend.justification = "top") +
    guides(fill = guide_legend(order = 1),
           color = guide_legend(order = 1),
           size = guide_legend(order = 2))
    
  ## Plot last year 
  ad1_coords <- 
    sf::st_point_on_surface(ad1_shape_simple) %>% 
    st_coordinates() %>% 
    data.table()
  ad1_coords$NAME <- ad1_shape_simple$admin1_name
  
  ad0_last_year <- 
    ad1_shape_simple %>% 
    left_join(admin1_est[year == max(year),], by =  c("admin1", "admin1_name")) %>% 
    mutate(sex_label = factor(sex, 1:2, c("Males", "Females"))) %>% 
    ggplot() + 
    geom_sf(aes(fill = mean)) +
    labs(fill = paste0("Mean\ndeaths","\n(", max(years), ")")) +
    geom_label_repel(data = ad1_coords,
                     size = 2.5,
                     aes(x = X, y = Y, label = NAME),
                     point.padding = unit(0.01, "lines"),
                     box.padding = unit(1.5, "lines"),
                     min.segment.length = unit(0, "lines"),
                     segment.alpha = 0.5) +
    scale_fill_distiller(palette = "RdYlBu",
                         direction = -1) + 
    labs(x = NULL, y = NULL, title = NULL) +
    coord_sf(datum = NA) +
    theme_classic(base_size = 10) + 
    facet_wrap(~ sex_label, nrow = 2)
  
  master_plot <- 
    arrangeGrob(ad1_time_trend, ad0_last_year,
                 layout_matrix = rbind(c(1, 1, 2),
                                       c(1, 1, 2)))
  
  pdf(file = paste0(model_dir, "/time_trend/all_admin1_plot.pdf"), height = 10, width = 15)
  grid.draw(master_plot)
  dev.off()
  
  # Make plot of age pattern of each Admin-1 ---------------------------------------------------------------------------------------------------
  message("Saving admin1 age graphs")
  ad1_code_list <- sort(unique(admin1_est_data$admin1))
  
  # Grab GBD data
  # Grab GBD estimates of HIV deaths
  source(paste0("<<<< FILEPATH REDACTED >>>>/lbd_hiv/sae/functions/vr_gbd_comparison_functions.r"))
  gbd_data <- load_gbd_estimates(ad1_code_list, years, cause_id = "298", include_data = T)
  
  pdf(file = paste0(model_dir, "/time_trend/admin1_age_trend_plots.pdf"), height = 10, width = 15)
  
  for (i in 1:length(ad1_code_list)) {
    adm1_code <- ad1_code_list[i]
    adm1_name <- unique(admin1_est_data[admin1 == adm1_code, admin1_name])
    pop_pct <- round(pop_admin1_pct[admin1 == adm1_code, percentile])
    
    # Subset estimates to admin1
    gbd_estimates <- data.table(gbd_data[["estimates"]])[loc_id == adm1_code]
    gbd_input     <- data.table(gbd_data[["data"]])[loc_id == adm1_code]
    
    adm1_data <- 
      admin1_data_complete %>% 
      filter(admin1 == adm1_code) %>% 
      left_join(gbd_input, by = c("year", "age", "sex")) %>% 
      rename(sae = events) %>% 
      mutate_at(vars(gbd), funs(.*100000)) %>% 
      gather(key = "Source", value = "data", sae, gbd) %>% 
      mutate(Source = factor(Source, c("sae", "gbd"), c("SAE", "GBD"))) %>% 
      mutate(sex_label = factor(sex, 1:2, c("Males", "Females")))
    
    gg_gbd <- 
      gbd_estimates %>% 
      left_join(admin1_est[admin1 == adm1_code],  c("age", "year", "sex")) %>% 
      mutate(sex_label = factor(sex, 1:2, c("Males", "Females"))) %>% 
      mutate_at(vars(gbd), funs(.*100000)) %>% 
      rename(sae = mean) %>% 
      gather(key = "Source", value = "mortality", sae, gbd) %>%
      mutate(Source = factor(Source, c("sae", "gbd"), c("SAE", "GBD"))) %>% 
      ggplot() + 
      geom_line(aes(x = year, y = mortality, color = sex_label, linetype = Source)) + 
      geom_point(data = adm1_data, aes(x = year, y = data, color = sex_label, shape = Source), alpha = 0.7) +
      scale_color_discrete(name = "Sex") +
      labs(y = "HIV deaths per 100,000 people", x = "Year",
           title = paste("HIV deaths in", adm1_name, "which has a population in", pop_pct, "percentile for", country_name)) +
      theme_light() + 
      facet_wrap(~ age, scales = "free")
    
    gg1 <- 
      admin1_est[admin1 == adm1_code] %>% 
      mutate(sex_label = factor(sex, 1:2, c("Males", "Females"))) %>% 
      ggplot() + 
      geom_point(data = filter(adm1_data, Source == "SAE"), aes(x = year, y = data, color = sex_label), shape = 1) +
      geom_line(aes(x = year, y = mean, color = sex_label)) + 
      geom_ribbon(aes(x = year, ymin = lb, ymax = ub, fill = sex_label), alpha = 0.2) + 
      scale_color_discrete(name = "Sex") +
      scale_fill_discrete(name = "Sex") +
      labs(title = paste("HIV deaths in", adm1_name, "which has a population in", pop_pct, "percentile for", country_name), 
           y = "HIV deaths per 100,000 people", x = "Year") +
      theme_light() + 
      facet_wrap(~ age, scales = "free")
    
    gg2 <- 
      admin1_est_data[admin1 == adm1_code] %>% 
      mutate(sex_label = factor(sex, 1:2, c("Males", "Females"))) %>% 
      ggplot() + 
      geom_point(aes(x = age, y = events, color = sex_label, size = pop), shape = 1) +
      geom_line(aes(x = age, y = mean, color = sex_label)) + 
      geom_ribbon(aes(x = age, ymin = lb, ymax = ub, fill = sex_label), alpha = 0.3) + 
      scale_color_discrete(name = "Sex") +
      scale_fill_discrete(name = "Sex") +
      labs(title = paste("HIV deaths in", adm1_name, "which has a population in", pop_pct, "percentile for", country_name), 
           y = "HIV deaths per 100,000 people", x = "Age group") +
      theme_light() + 
      facet_wrap(~ year)
    
    # source(paste0(core_repo, "/mbg_central/visualization_functions.R"))
    # admin1_loc <- location_map_draw(as(filter(ad1_shape_simple, admin1 == adm1_code), "Spatial"), as(ad0_shape_simple, "Spatial"))
    # 
    
    # plot1 <- 
    #   arrangeGrob(gg1, admin1_loc,
    #               layout_matrix = rbind(c(1, 1, 1, 1, 2),
    #                                     c(1, 1, 1, 1, 2),
    #                                     c(1, 1, 1, 1, 2)))
    # 
    # plot2 <- 
    #   arrangeGrob(gg2, admin1_loc,
    #               layout_matrix = rbind(c(1, 1, 1, 1, 2),
    #                                     c(1, 1, 1, 1, 2),
    #                                     c(1, 1, 1, 1, 2)))
    print(gg_gbd)
    # grid.draw(plot1)
    # plot.new()
    # grid.draw(plot2)
    #if (i != length(ad1_code_list)) plot.new()
  }
  dev.off()
}


##############################################
# Model results compare for SAE models
# This function plots two model results to compare model runs
##############################################

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(dplyr)
library(purrr)
library(tidyr)
library(sf)

model_results_compare_sae <- function(country, years, save_file, covariate_data = NULL, 
                                      run_dates = NULL, run_date_labels = NULL) {
  
  library(gridExtra)
  library(grid)
  select <- dplyr::select
  rename <- dplyr::rename
  
  if (is.null(covariate_data)) {
    setwd("<<<< FILEPATH REDACTED >>>>/lbd_core/sae_central/")
    source("settings.r")
    
    ## Load final model output data and shapefiles -----------------------------------------------------
    # Estimates
    data <- 
      Reduce(merge, (lapply(1:length(run_dates), function(i) {
        main_dir <-  paste0("<<<< FILEPATH REDACTED >>>>")
        load(paste0(main_dir,"/est_all.rdata"))
        est <- est[sex == 3 & age == 99,]
        est <- est[, .(area, year, mean)]
        
        # Convert into rate per 100,000
        est[, mean := mean*100000]
        setnames(est, "mean", run_date_labels[i])
      })))
  } else {
    run_dates <- as.character(years)
    run_date_labels <- as.character(years)
    if (length(names(covariate_data)) != 3 | !"area" %in% names(covariate_data) | !"year" %in% names(covariate_data)){
      stop("Need to include covariate data frame with only columns of area, year and then the covariate measure")
    }
    cov_name <- setdiff(names(covariate_data), c("year", "area"))
    data <- 
      covariate_data %>% 
      spread(year, get(cov_name)) %>% 
      data.table()
  }
  
  for (i in 1:length(run_dates)) {
    for (j in 1:length(run_dates)) {
      if (i == j) next
      data[, paste(run_date_labels[i], run_date_labels[j]) := get(paste(run_date_labels[i])) - get(paste(run_date_labels[j]))]
    }
  }
  
  if (is.null(covariate_data)) data_long <- melt(data, id.vars = c("area", "year")) else data_long <- melt(data, id.vars = c("area"))
  data_long[, type := ifelse(variable %in% run_date_labels, 'level', 'comp')]
  
  message("Making maps & plots")
  # set up a list of all maps & scatters to make, along with labels
  all_plots <- CJ(mod1 = run_date_labels, mod2 = run_date_labels, sorted = F)
  all_plots[mod1 == mod2, `:=` (label = mod1, variable = mod1)]
  all_plots[mod1 != mod2, `:=` (label = paste(mod1, '-', mod2), variable = paste(mod1, mod2))]
  all_plots[, type := ifelse(mod1 == mod2, "level", "comp")]
  
  # Makes scatterplot in lower 
  all_plots[, geom := ifelse(as.numeric(factor(mod1, unique(mod1))) < as.numeric(factor(mod2, unique(mod2))), 'scatter', 'map')]
  # get limits and color scales across all plots
  range <- data_long[, list(min = min(value), max = max(value)), keyby = type]
  if (range['level', min] < 0) range['level', c('min', 'max') := as.list(c(-1, 1) * max(abs(c(min, max))))]
  
  brks <- list(level = pretty(range['level', c(min, max)], n = 5),
               comp = pretty(range['comp', c(min, max)], n = 5))
  
  colors <- list(level = if (range['level', min] < 0) brewer.pal(7, 'PiYG') else brewer.pal(8, 'RdPu')[-1],
                 comp = brewer.pal(7, 'PuOr'))
  
  values <- data_long[, as.list(quantile(value, c(0.01, 0.99))), keyby = type]
  if (values['level', `1%`] < 0) values['level', c('1%', '99%') := as.list(c(-1, 1) * max(abs(c(`1%`, `99%`))))]
  values <- 
    sapply(c('level', 'comp'), function(x) {
      rescale(c(range[x, min], seq(values[x, `1%`], values[x, `99%`], length.out = 5), range[x, max]))
    }, simplify = F)
  
  # make all maps in all years
  plot_years <- sapply(years, paste, collapse = "-")
  if (!is.null(covariate_data)) plot_years = 1
  plot_list <- sapply(plot_years, function(y) {
    lapply(1:nrow(all_plots), function(p) {
      this_geom <- all_plots[p, geom]
      this_type <- all_plots[p, type]
      
      if (this_geom == 'map') {
        # Grab shapefile and append correct data
        shape <- 
          st_read(vr_pull_shp(paste0(country, "_", admin), "stable", "full")$shapefile, quiet = T) %>% 
          mutate(uid = as.numeric(as.character(uid)))
        if (is.null(covariate_data)) fdata = data[year == y,] else fdata = data
        fdata = fdata[, .(area, get(all_plots[p, variable]))]
        
        shape %>% 
          left_join(fdata, by = c("uid" = "area")) %>% 
          rename(value = V2) %>% 
          ggplot() + 
          geom_sf(aes(fill = value), size = .3, color = NA) + 
          scale_fill_gradientn(colors = colors[[this_type]], values = values[[this_type]],
                               limits = range(brks[[this_type]]), breaks = brks[[this_type]],
                               name = if (this_type == "level") "Level" else "Diff.") + 
          labs(title = all_plots[p, label]) +
          theme_void() + 
          coord_sf(datum = NA) +
          theme(plot.title = element_text(hjust=0.5), plot.margin = unit(c(0, 0, 0, 0), "cm")) +
          guides(fill = guide_colorbar(barwidth = 0.5, barheight = 8))
        
        
      } else if (this_geom == 'scatter') {
        cuts <- seq(min(brks$level), max(brks$level), length.out = 100)
        mids <- (cuts[-1] + cuts[-100])/2
        if (is.null(covariate_data)) fdata = data[year == y,] else fdata = data
        fdata <- data[, list(x = get(all_plots[p, mod2]), y = get(all_plots[p, mod1]))]
        fdata[, x := mids[cut(x, cuts, labels=F)]]
        fdata[, y := mids[cut(y, cuts, labels=F)]]
        fdata <- fdata[, .N, keyby='x,y']
        
        ggplot(fdata, aes(x = x, y = y, size = N)) +
          geom_point(alpha = 0.7, show.legend = F) +
          geom_abline(intercept = 0, slope = 1) +
          scale_size_area() +
          coord_equal(ratio = 1, xlim = range(brks$level), ylim = range(brks$level)) +
          labs(x = all_plots[p, mod2], y = all_plots[p, mod1]) +
          theme_classic() +
          theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
        
      }
    })
  }, simplify = F)
  
  ## Plot and save maps ----------------------------------------------------------------------------
  message("Saving")
  
  
  if (sum(run_dates != run_date_labels) > 0) {
    cap <- paste("Models:", paste(paste0(run_dates, " (", run_date_labels, ")"), collapse = ", "))
  } else {
    cap <- paste("Models:", paste(run_dates, collapse = ", "))
  }
  
  pdf(save_file, width = 4*length(run_dates), height = 4*length(run_dates) + 1.5)
  
  if (is.null(covariate_data)){
    for (y in plot_years) {
      grid.arrange(grobs = plot_list[[as.character(y)]],
                   top = textGrob(y, gp=gpar(cex=2), just="top"),
                   bottom = textGrob(cap, gp=gpar(cex=1), just="bottom"))
    }
  } else {
    grid.arrange(grobs = plot_list[[1]],
                 top = textGrob(cov_name, gp=gpar(cex=2), just="top"),
                 bottom = textGrob("", gp=gpar(cex=1), just="bottom"))
  }
  
  dev.off()
  
  return(paste("Plots saved:", save_file))
}


