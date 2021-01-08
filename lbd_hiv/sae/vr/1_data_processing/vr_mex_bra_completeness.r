library(tidyr)
library(stringr)
library(scales)
library(foreign)
library(MCMCglmm)
library(readr)
library(parallel)
library(fitdistrplus)
library(tidybayes, lib.loc = "<<<< FILEPATH REDACTED >>>>")

rm(list = ls())

## Print message for clarity
message(commandArgs()[4])
message(commandArgs()[5])
message(commandArgs()[6])
message(commandArgs()[7])

## Set indicator
country         <- commandArgs()[4]
rerun           <- commandArgs()[5]
calc_child_comp <- commandArgs()[6]
calc_adult_comp <- commandArgs()[7]

# Setup -------------------------------------------------------------------------------------------------------------
# Load LBD functions needed
core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>")
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv")) %>% setdiff(("plyr"))
mbg_setup(package_list = package_list, repos = core_repo)
source(paste0(core_repo, "/sae_central/settings.r"))
source(paste0(core_repo, "/data_central/vr_functions.R"))
source(paste0(core_repo, "/data_central/sae_functions.R"))

# Load HIV specific functions
function_files = list.files(paste0("<<<< FILEPATH REDACTED >>>>/lbd_hiv/sae/functions/"), full.names = T)
sapply(function_files, source)

# Source central functions
source("<<<< FILEPATH REDACTED >>>>/r/get_age_metadata.R")
source("<<<< FILEPATH REDACTED >>>>/r/get_population.R")

# Avoid masking dplur, setwd
select    <- dplyr::select
summarize <- dplyr::summarize
setwd("<<<< FILEPATH REDACTED >>>>/lbd_core/sae_central/")

# Assign number of draws
n_draws = 1000
gbd_locs <- get_gbd_locs(country, shapefile_version = "current", rake_subnational = T)

# Grab most recent data input settings
data_run_date <- most_recent_date(paste0("<<<< FILEPATH REDACTED >>>>"))
settings_folder <- paste0("<<<< FILEPATH REDACTED >>>>")
get_settings(settings_folder)

## Adjusting geographic variation to national VR completeness for adults ------------------------------------------------------
if (calc_adult_comp) {
  
  # Grab national draws of adult completeness from DDM methods
  library(rhdf5)
  subnat_completeness_draws <-
    rbindlist(mclapply(0:{n_draws - 1}, mc.cores = 20, function(draw) {
      draw = as.character(draw)
      data <- 
        h5read("<<<< FILEPATH REDACTED >>>>", draw) %>% 
        filter(year %in% years, location_id %in% gbd_locs$location_id) %>%
        mutate(completeness = 10 ^ (trunc_pred)) %>%
        dplyr::select(location_id, year, draw, completeness) %>%
        data.table()
    }))
  
  
  # Have draws ordered from 1 to 1000
  subnat_completeness_draws[, draw := draw + 1]
  
  sample_rows <- nrow(subnat_completeness_draws[completeness == 1])
  subnat_completeness_draws[completeness == 1, completeness := runif(sample_rows, min = 0.95, max = 1)]
  
  # Make final completeness
  final_adult_completeness <- copy(subnat_completeness_draws)
  setnames(final_adult_completeness, c("location_id", "completeness"), c("area", "final_complete"))
  saveRDS(final_adult_completeness, file = paste0("<<<< FILEPATH REDACTED >>>>"))
  
  
  areas = sort(unique(final_adult_completeness$area))
  years = sort(unique(final_adult_completeness$year))
  
  adult_prior_dist <-
    rbindlist(mclapply(areas, mc.cores = 10, function(a) {
      message("working on area ", a)
      rbindlist(lapply(years, function(y) {
        logit_norm_fit <-
          final_adult_completeness %>%
          filter(area == a, year == y) %>%
          mutate(final_complete = logit(final_complete)) %>%
          pull(final_complete) %>%
          fitdist(dist = "norm", method = "mme")
        
        fit_table <- data.table(area = a,
                                year = y,
                                mean = logit_norm_fit$estimate[["mean"]],
                                sd = logit_norm_fit$estimate[["sd"]])
      }))
    }))
  
  # Make graphs of prior fit 
  gg_graphs_logit <- 
    lapply(areas, function(a) {
      message("working on area ", a)
    
    # Simulate from prior districution
      area_prior <-
        adult_prior_dist[area == a, .(normal_prior = inv.logit(rnorm(n_draws, mean = mean, sd = sd)),
                                      draw = 1:n_draws), by = c("area", "year")]
      
      # Merge on area specific completeness
      area_completeness <- merge(final_adult_completeness[area == a], area_prior, by = c("area", "year", "draw"))
    
      int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0]
      gg_complete <-
        area_completeness %>%
          ggplot() +
          stat_halfeye(aes(x = year, y = final_complete), side = "left") +
          stat_slab(aes(x = year, y = normal_prior), side = "left", color = "red", fill = NA, linetype = "dotted") +
          labs(y = "completeness", x = "location_id") +
          scale_x_continuous(breaks = int_breaks) +
          theme_classic() +
          coord_cartesian(ylim = c(0, 1)) +
          facet_wrap(~ area)
    })
  
  total_areas <- length(areas)
  indices <- data.table(value1 = seq(1, total_areas, 4), value2 = c(seq(4, total_areas, 4), total_areas))
  
  # Make graph of priors, batched by 10 areas to one sheet
  save_file = paste0(settings_folder, "<<<< FILEPATH REDACTED >>>>")
  pdf(file = save_file, height = 10, width = 18)
  apply(indices, 1, function(ind){
    rows = ceiling(length(ind[[1]]:ind[[2]]) / 2) 
    do.call(grid.arrange, gg_graphs_logit[ind[[1]]:ind[[2]]])
  })
  dev.off() 
  saveRDS(adult_prior_dist, file = paste0("<<<< FILEPATH REDACTED >>>>"))
  
} else {
  adult_prior_dist <- data.table(expand.grid(area = unique(gbd_locs$location_id), year = years))[, mean := 10][, sd := 0.01]
}

if (calc_child_comp) {
  
  # Draws of under 5 deaths from GBD 2019 -- May need to update
  source("<<<< FILEPATH REDACTED >>>>/r/get_draws.R")
  gbd_draws <-
    get_draws('cause_id',
              measure_id = 1,
              source = 'codcorrect',
              metric_id = 1,
              gbd_id = 294,
              year_id = years,
              gbd_round_id = 5,
              #version_id = 135,
              location_id = gbd_locs$location_id,
              age_group_id = 39, # 0-14 year olds
              sex_id = 3)
              #decomp_step = 'step4')
  
  # Sum over all years 
  cols <- paste0("draw_", seq(0, n_draws-1, 1))  # names of columns 
  draws_child_deaths <- 
    gbd_draws %>% 
    pivot_longer(cols = cols, names_to = "draw") %>% 
    mutate(draw = as.numeric(str_remove(draw, "draw_"))) %>% 
    mutate(draw = draw + 1) %>% # draws from 1 to 1000
    dplyr::select(area = location_id, year = year_id, age = age_group_id, draw, value) %>% 
    group_by(area, year, draw) %>% 
    dplyr::summarize(deaths = sum(value)) %>% 
    ungroup() %>% 
    arrange(area, year, draw) %>% 
    data.table() 
  
  # Calculate child deaths at area and national level
  vr_deaths <- vr_prep_all_cause_mortality(country, "adm2", years, agg_level= "admin1")
  area_child_deaths <- vr_deaths[age < 15, .(deaths = sum(deaths)), by = .(area, year)]
  setnames(area_child_deaths, "deaths", "reg_deaths")
  
  # Get completeness draws for admin1
  draws_child_comp <- 
    draws_child_deaths %>% 
    left_join(area_child_deaths, by = c("area", "year")) %>% 
    mutate(completeness = reg_deaths / deaths) %>% 
    data.table()
  
  int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0]
  gg_complete <-
    draws_child_comp %>%
    ggplot() +
    stat_halfeye(aes(x = year, y = completeness), side = "left") +
    labs(y = "child completeness", x = "location_id") +
    scale_x_continuous(breaks = int_breaks) +
    theme_classic() +
    facet_wrap(~ area)
  
  pdf(file = paste0(settings_folder, "<<<< FILEPATH REDACTED >>>>"), height = 10, width = 18)
  print(gg_complete)
  dev.off()
  
  # See how many draws of child completeness are above 1 
  message(round(100 * nrow(draws_child_comp[completeness > 1]) / nrow(draws_child_comp), 2), 
          "% of child completeness draws are above 1")
  
  
  # Make sure child completeness is below 1, or nudge registered deaths
  if (nrow(draws_child_comp[completeness > 1]) != 0) {
    # Make registered deaths a little less than deaths using same completeness cutoff
    draws_child_comp[completeness > 1, deaths := reg_deaths/(0.99)]
    draws_child_comp[, completeness := reg_deaths / deaths]
    draws_child_deaths <- draws_child_comp[, .(year, draw, deaths)]
  }
  
  # Assert that completeness is now less than 1
  assert_that(
    draws_child_comp %>%
      mutate(completeness = reg_deaths / deaths) %>%
      filter(completeness >= 1) %>%
      nrow() == 0
  )
  
  # Make final completeness
  final_child_completeness <- copy(draws_child_comp)
  setnames(final_child_completeness, "completeness", "final_complete")
  saveRDS(final_child_completeness, file = paste0("<<<< FILEPATH REDACTED >>>>"))
  
  areas = sort(unique(final_child_completeness$area))
  years = sort(unique(final_child_completeness$year))
  
  # make final completeness
  child_prior_dist <-
    rbindlist(mclapply(areas, mc.cores = 5, function(a) {
      message("working on area ", a)
      rbindlist(lapply(years, function(y) {
        logit_norm_fit <-
          final_child_completeness %>%
          filter(area == a, year == y) %>%
          mutate(final_complete = logit(final_complete)) %>%
          pull(final_complete) %>%
          fitdist(dist = "norm", method = "mme")
        
        fit_table <- data.table(area = a,
                                year = y,
                                mean = logit_norm_fit$estimate[["mean"]],
                                sd = logit_norm_fit$estimate[["sd"]])
      }))
    }))
  
  # Make graphs of prior fit 
  gg_graphs_logit <- 
    lapply(areas, function(a) {
      message("working on area ", a)
      
      # Simulate from prior districution
      area_prior <-
        child_prior_dist[area == a, .(normal_prior = inv.logit(rnorm(n_draws, mean = mean, sd = sd)),
                                      draw = 1:n_draws), by = c("area", "year")]
      
      # Merge on area specific completeness
      area_completeness <- merge(final_child_completeness[area == a], area_prior, by = c("area", "year", "draw"))
      
      int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0]
      gg_complete <-
        area_completeness %>%
        ggplot() +
        stat_halfeye(aes(x = year, y = final_complete), side = "left") +
        stat_slab(aes(x = year, y = normal_prior), side = "left", color = "red", fill = NA, linetype = "dotted") +
        labs(y = "completeness", x = "location_id") +
        scale_x_continuous(breaks = int_breaks) +
        theme_classic() +
        coord_cartesian(ylim = c(0, 1)) +
        facet_wrap(~ area)
    })
  
  total_areas <- length(areas)
  indices <- data.table(value1 = seq(1, total_areas, 4), value2 = c(seq(4, total_areas, 4), total_areas))
  
  # Make graph of priors, batched by 10 areas to one sheet
  save_file = paste0(settings_folder, "/figures/child_completeness_prior.pdf")
  pdf(file = save_file, height = 10, width = 18)
  apply(indices, 1, function(ind){
    rows = ceiling(length(ind[[1]]:ind[[2]]) / 2) 
    do.call(grid.arrange, gg_graphs_logit[ind[[1]]:ind[[2]]])
  })
  dev.off() 
  
  saveRDS(child_prior_dist, file = paste0("<<<< FILEPATH REDACTED >>>>"))
  
} else {
  child_prior_dist <- data.table(expand.grid(area = unique(gbd_locs$location_id), year = years))[, mean := 10][, sd := 0.01]
}

completeness_prior_dist <- 
  rbind(adult_prior_dist[, age_group := 1], child_prior_dist[, age_group := 0])[order(age_group, area, year)]

saveRDS(completeness_prior_dist, file = paste0(settings_folder, "/comp_prior_dist.RDS"))



