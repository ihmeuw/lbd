# Script to generate completeness priors for Colombia, Ecuador and Guatemala, Mexico or Brazil 

library(tidyr)
library(stringr)
library(scales)
library(foreign)
library(MCMCglmm)
library(readr)
library(parallel)
library(fitdistrplus)

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
package_list <- readLines("<<<< FILEPATH REDACTED >>>>/package_list.csv") %>% setdiff(("plyr"))
mbg_setup(package_list = package_list, repos = core_repo)
source(paste0(core_repo, "/sae_central/settings.r"))
source(paste0(core_repo, "/data_central/vr_functions.R"))
source(paste0(core_repo, "/data_central/sae_functions.R"))

# Load HIV specific functions
function_files = list.files(paste0("<<<< FILEPATH REDACTED >>>>/lbd_hiv/sae/functions/"), full.names = T)
sapply(function_files, source)

# Source central functions
source("<<<< FILEPATH REDACTED >>>>/get_age_metadata.R")
source("<<<< FILEPATH REDACTED >>>>/get_population.R")

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

# Grab shapefile
shape <-
  st_read(vr_pull_shp(paste0(country, "_adm2"), "stable", "full")$shapefile, quiet = T) %>%
  mutate(uid = as.numeric(as.character(uid))) %>%
  arrange(uid) %>% select(uid)
shapefile_path <- vr_pull_shp(paste0(country, "_adm2"), "stable", "full")$shapefile


if (rerun) {
  
  # Check completeness estimates from DDM if needed
  gbd_locs <- get_gbd_locs(country, shapefile_version = "current", rake_subnational = T)
  
  ## Underlying geographic variation in VR completeness ----------------------------------------------------------------------
  # Get numerator of initial completeness, which is the smoothed number of all-cause deaths in each area across all years 
  vr_deaths_list <- completeness_shrink_u5m_deaths(country, years, pop_file, shape, adjmat_file)
  shrinkage_model <- vr_deaths_list[[1]]
  adjusted_vr <- vr_deaths_list[[2]]
  
  # Save vr deaths to scratch
  dir.create(paste0("<<<< FILEPATH REDACTED >>>>"), recursive = T)
  saveRDS(vr_deaths_list, file = paste0("<<<< FILEPATH REDACTED >>>>"))
  
  # Work on denominator of completeness, load draws of U5M deaths per pixel and aggregate
  if (country %in% c("gtm", "cri")) reg <- "caca-mex" else reg <- "ansa+trsa-bra"
  modeling_shapefile_version <- "2019_05_06"
  
  # Get draws of the denominator, which are the number of deaths in those same areas
  death_draws <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
  
  # Make simple raster for region
  subset_shape <-
    load_simple_polygon(gaul_list = get_adm0_codes(reg, shapefile_version = modeling_shapefile_version),
                        buffer = 0.4,
                        subset_only = FALSE,
                        shapefile_version = modeling_shapefile_version)[[1]]
  reg_raster <- build_simple_raster_pop(subset_shape)[["simple_raster"]]
  
  # Set the number of draws you want to use, to then be fractionally aggregated
  n_draws <- 1000
  cols <- paste0("draw_", seq(1, n_draws, 1))
  
  # Make raster of each draw from cell pred
  yrs = dim(death_draws)[1]/length(cellIdx(reg_raster))
  message(paste("Making a RasterBricks with", yrs, "layers"))
  
  draw_rasters <- 
    lapply(seq(1, n_draws, 1), function(draw){
      message("Making rasters for draw ", draw)
      draw_raster <- insertRaster(reg_raster,  matrix(death_draws[, draw],  ncol = yrs))
      draw_raster <- draw_raster[[years - 1999]] # subset to correct years
      return(draw_raster)
    })
  names(draw_rasters) <- cols
  
  # Get fractional aggregation for each draw
  cov_config <- data.table(covariate = cols, agg_method = "sum", release = NA)
  lbd_draws <- 
    frac_agg_covs(cov_config = cov_config,
                  years = years,
                  shapefile_path = shapefile_path,
                  shapefile_field = "uid",
                  core_repo = core_repo,
                  cores = 10, 
                  custom_covariate = draw_rasters)
  setnames(lbd_draws, "uid", "area")
  
  # Make draws longer, combine deaths across all years in each area
  lbd_draws <- 
    lbd_draws %>% 
    pivot_longer(cols = cols, names_to = "draw") %>% 
    mutate(draw = as.numeric(str_remove(draw, "draw_"))) %>% 
    group_by(draw, area) %>% 
    summarize(value = sum(value)) %>% 
    ungroup() %>% 
    data.table()
  
  # Save draws for to sctach
  saveRDS(lbd_draws, file = paste0("<<<< FILEPATH REDACTED >>>>"))
  
  ## Calculate initial completeness distribution -------------------------------------------------------
  # Initial completeness numerator
  vr_deaths <- 
    readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))[["adjusted_vr"]] %>% 
    dplyr::select(area, vr = adjusted)
  
  # Initial completeness denominator
  lbd_draws <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
  if ("draws" %in% names(lbd_draws)) setnames(lbd_draws, "draws", "draw") # for consistency
  
  # For colombia add missing islands from lbd draws, TEMP FIX
  if (country == "col") {
    # Grab draws of mortality rate
    died <- fread("<<<< FILEPATH REDACTED >>>>")[, ADM2_CODE := name][, -1]
    
    # Grab conversion rate 
    conversion <- fread("<<<< FILEPATH REDACTED >>>>")
    
    # Pop an NAX are used to convert to deaths
    load(pop_file)
    pop <- pop[age == 0, .(pop = sum(pop)), by = .(area, year)]
    
    # Draws of deaths
    islands <- 
      died %>% 
      filter(ADM2_CODE %in% c(3026050, 2026050)) %>% 
      pivot_longer(cols = paste0("V", 1:250), names_to = "draws") %>% 
      mutate(draw = as.numeric(str_remove(draws, "V"))) %>% 
      mutate(area = ifelse(ADM2_CODE == 3026050, 1076, 1077)) %>% # San Andrew is 3026050 and 1076
      select(area, year, draw, value) %>% 
      left_join(conversion[, .(year = year_id, a)], by = "year") %>% 
      left_join(pop, by = c("area", "year")) %>% 
      mutate(value = pop * (1 / (5 / value + a - 5))) %>% # Have deaths by year
      group_by(area, draw) %>% # Now summarize over all years for completeness estimate
      summarize(value = sum(value)) %>% 
      ungroup() %>% 
      data.table()
    
    # Repeat each row 4 times to get to 1000 draws
    island_draws <- 
      islands %>%
      mutate(freq = 4) %>%
      uncount(freq) %>%
      group_by(area) %>% 
      mutate(draw = 1:1000) %>% 
      ungroup()
      data.table()
    
    lbd_draws <- rbind(lbd_draws, island_draws)
  }
  
  # Draws of initial completeness, vr deaths over the draws of u5m deaths
  initial_completeness <- merge(lbd_draws, vr_deaths, by = "area")[, initial_complete := vr / value][order(area, draw)]
  
  # For now, cap completeness at 1. These are now geographic distriution of completeness
  message("Percent of draws above 1 is ",
          round(100*nrow(initial_completeness[initial_complete >= 1]) / nrow(initial_completeness), 1), "%")
  
  # Truncate initial completeness
  # Current truncate to 99 percentile of value under 1 or 0.99
  initial_completeness[, initial_complete := ifelse(initial_complete >= 1, 0, initial_complete)]
  initial_completeness[, initial_complete := ifelse(initial_complete == 0, max(initial_complete), initial_complete), by = .(area)]
  initial_completeness[, initial_complete := ifelse(initial_complete == 0, max(c(quantile(initial_complete, 0.99), 0.99)), initial_complete), by = .(area)]
  
  # If all estimates were above 1, replace with maximum value
  initial_completeness[, initial_complete := ifelse(initial_complete == 0, max(initial_complete), initial_complete)]
  initial_completeness <- initial_completeness[, .(draw, area, initial_complete)]
  
  ## Adjusting geographic variation to national VR completeness for adults ------------------------------------------------------
  if (calc_adult_comp) {
    
    # Grab national draws of adult completeness from DDM methods
    library(rhdf5)
    nat_completeness_draws <-
      rbindlist(mclapply(0:{n_draws - 1}, mc.cores = 25, function(draw) {
        draw = as.character(draw)
        data <- 
          h5read("<<<< FILEPATH REDACTED >>>>", draw) %>% 
          filter(year %in% years, location_id %in% gbd_locs$location_id) %>%
          mutate(completeness = 10 ^ (trunc_pred)) %>%
          dplyr::select(location_id, year, draw, completeness) %>%
          data.table()
      }))
    
    
    # Have draws ordered from 1 to 1000
    nat_completeness_draws[, draw := draw + 1]
    
    # Some completeness draws have 5 versions of each draw, for those take a random one of those draws
    # This is because of the 5 DDM methods in certain years 
    if (nrow(nat_completeness_draws[, .(number = .N), by = .(year, draw)][number > 1]) != 0) {
      five_draw_year <- nat_completeness_draws[, .(number = .N), by = .(year, draw)][number > 1, unique(year)]
      five_draw_completeness <- nat_completeness_draws[year == five_draw_year, .(completeness= nth(completeness, sample(1:5, 1000, replace = T)[draw])), .(location_id, year, draw)]
      nat_completeness_draws <- rbind(five_draw_completeness, nat_completeness_draws[year != five_draw_year])[order(year, draw)]
    }
    
    # Distribute completeness that is exactly 1 to between 0.95 and 1 given that is the procedure taken by GBD
    sample_rows <- nrow(nat_completeness_draws[completeness == 1])
    nat_completeness_draws[completeness == 1, completeness := runif(sample_rows, min = 0.95, max = 1)]
    
    # Plot distributions of national completeness
    gg <- 
      ggplot(nat_completeness_draws) + 
      geom_histogram(aes(x = completeness)) + 
      theme_bw() + 
      labs(x = "draws of national completeness", title = paste("Completeness draws for", country)) + 
      facet_wrap(~ year)
    ggsave(filename = paste0(settings_folder, "/figures/ddm_draws.pdf"), gg)
    saveRDS(nat_completeness_draws, file = paste0("<<<< FILEPATH REDACTED >>>>"))
    
    # Load draws of adult completeness 
    nat_completeness_draws <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
    
    # Calculate adult deaths at area and national level
    vr_deaths <- vr_prep_all_cause_mortality(country, "adm2", years)
    area_adult_deaths <- vr_deaths[age >= 15, .(deaths = sum(deaths)), by = .(area, year)]
    
    # Draws of national deaths uses draws of national completeness
    national_adult_deaths <- vr_deaths[age >= 15, .(deaths = sum(deaths)), by = .(year)]
    draws_national_adult_deaths <- 
      merge(nat_completeness_draws, national_adult_deaths, by = "year") %>% 
      mutate(deaths = deaths / completeness) %>% 
      select(year, draw, deaths) %>% 
      data.table()
    
    # Rake such that area level completeness is equal to national level
    message("Finding logit raking factors adult completeness")
    adult_logit_rf <-
      rbindlist(lapply(years, function(y) {
        message(paste("Working on year", y, "->>>>"))
        rbindlist(mclapply(1:1000, mc.cores = 50, function(d) {
          
          # Grab GBD total deaths
          gbdval <-  draws_national_adult_deaths[year == y & draw == d, deaths]
          
          # Draw of area level completeness
          pixelval  <- matrix(initial_completeness[draw == d, initial_complete], ncol = 1)
          
          # Weighted by number of VR deaths per area in that year
          weightval <- matrix(area_adult_deaths[year == y, deaths], ncol = 1)
          
          r <- 
            LogitFindK_mort(gbdval     = gbdval,
                            pixelval   = pixelval,
                            weightval  = weightval,
                            approx_0_1 = T,
                            MaxIter    = 100)
          
          # Return as data table
          return(data.table(year = y,
                            rf = r,
                            draw = d))
        }))
      }))
    
    # Apply correction to get final completeness 
    final_adult_completeness <- 
      initial_completeness %>% 
      left_join(adult_logit_rf, by = "draw") %>% 
      mutate(final_complete = invlogit(logit(initial_complete) + rf)) %>% 
      data.table()
    
    # Save final completeness -- this contains the adult rf 
    saveRDS(final_adult_completeness, file = paste0("<<<< FILEPATH REDACTED >>>>"))
    
    # Fit distribution to logitnormal priors pased in as arguments
    library(fitdistrplus)
    adult_prior_dist <-
      rbindlist(mclapply(sort(unique(final_adult_completeness$area)), mc.cores = 10, function(a) {
        message("working on area ", a)
        rbindlist(lapply(sort(unique(final_adult_completeness$year)), function(y) {
          logit_norm_fit <-
            final_adult_completeness %>%
            filter(area == a, year == y) %>%
            mutate(final_complete = logit(final_complete)) %>%
            pull(final_complete) %>%
            fitdist(dist = "norm", method = "mle")
          
          fit_table <- data.table(area = a,
                                  year = y,
                                  mean = logit_norm_fit$estimate[["mean"]],
                                  sd = logit_norm_fit$estimate[["sd"]])
        }))
      }))
    
    saveRDS(adult_prior_dist, file = paste0("<<<< FILEPATH REDACTED >>>>"))
    
    # Make visualization plots, saved in figures folder 
    dir.create(paste0(settings_folder, "/figures/"))
    plot_completeness_prior(final_adult_completeness, adult_prior_dist, paste0(settings_folder, "/figures/adult_prior_dist.pdf"))
    plot_dist_comparison(final_adult_completeness, adult_prior_dist, paste0(settings_folder, "/figures/adult_prior_compare.pdf"))
    
  } else {
    # Make dummy mean and SD since this will not be used
    adult_prior_dist <- data.table(expand.grid(area = unique(initial_completeness$area), year = years))[, mean := 10][, sd := 0.01]
  }
  
  ## Now calculate child completeness -------------------------------------------------------------------------------------
   if (calc_child_comp) {
  
    # Draws of under 5 deaths from GBD 2019 -- May need to update
    source("<<<< FILEPATH REDACTED >>>>/get_draws.R")
    gbd_draws <-
      get_draws('cause_id',
                measure_id = 1,
                source = 'codcorrect',
                metric_id = 1,
                gbd_id = 294,
                year_id = years,
                gbd_round_id = 5,
                location_id = gbd_locs$location_id,
                age_group_id = 39, # 0-14 year olds
                sex_id = 3)
    
    
      # get_draws('cause_id',
      #           measure_id = 1,
      #           source = 'codcorrect',
      #           metric_id = 1,
      #           gbd_id = 294,
      #           year_id = years,
      #           gbd_round_id = 6,
      #           version_id = 135,
      #           location_id = gbd_locs$location_id,
      #           age_group_id = 39, # 0-14 year olds
      #           sex_id = 3,
      #           decomp_step = 'step4')
    
    # Sum over all years 
    cols <- paste0("draw_", seq(0, n_draws-1, 1))  # names of columns 
    draws_national_child_deaths <- 
      gbd_draws %>% 
      pivot_longer(cols = cols, names_to = "draw") %>% 
      mutate(draw = as.numeric(str_remove(draw, "draw_"))) %>% 
      mutate(draw = draw + 1) %>% # draws from 1 to 1000
      dplyr::select(year = year_id, age = age_group_id, draw, value) %>% 
      group_by(year, draw) %>% 
      dplyr::summarize(deaths = sum(value)) %>% 
      ungroup() %>% 
      arrange(year, draw) %>% 
      data.table() 
    
    # Calculate child deaths at area and national level
    vr_deaths <- vr_prep_all_cause_mortality(country, "adm2", years)
    area_child_deaths <- vr_deaths[age < 15, .(deaths = sum(deaths)), by = .(area, year)]
    
    # Plot draws of national child completeness to make sure none are above 1
    national_child_deaths <- vr_deaths[age < 15, .(reg_deaths = sum(deaths)), by = .(year)]
    
    # Get completeness draws
    draws_national_child_comp <- 
      draws_national_child_deaths %>% 
      left_join(national_child_deaths, by = "year") %>% 
      mutate(completeness = reg_deaths / deaths) %>% 
      data.table()
    
    gg_child_comp <- 
      draws_national_child_comp %>% 
      ggplot() + 
      geom_histogram(aes(x = completeness)) + 
      geom_vline(xintercept = 1) + 
      labs(x = "draws of completeness") + 
      theme_bw() + 
      facet_wrap(~ year) 
    
    pdf(file = paste0(settings_folder, "/national_child_completeness.pdf"), height = 10, width = 18)
    gg_child_comp
    dev.off()
    
    # Make sure child completeness is below 1, or nudge registered deaths
    if (nrow(draws_national_child_comp[completeness > 1]) != 0) {
      # Make registered deaths a little less than deaths using same completeness cutoff
      draws_national_child_comp[completeness > 1, deaths := reg_deaths/(0.99)]
      draws_national_child_comp[, completeness := reg_deaths / deaths]
      draws_national_child_deaths <- draws_national_child_comp[, .(year, draw, deaths)]
    }
    
    # Assert that completeness is now less than 1
    assert_that(draws_national_child_deaths %>% 
                  left_join(national_child_deaths, by = "year") %>% 
                  mutate(completeness = reg_deaths / deaths) %>%
                  filter(completeness >= 1) %>% 
                  nrow() == 0)
    
    
    # Rake such that area level completeness is equal to national level
    message("Finding logit raking factors for child completeness")
    child_logit_rf <-
      rbindlist(lapply(years, function(y) {
        message(paste("Working on year", y, "->>>>"))
        rbindlist(mclapply(1:1000, mc.cores = 50, function(d) {
          
          # Grab GBD total deaths
          gbdval <-  draws_national_child_deaths[year == y & draw == d, deaths]
          
          # Draw of area level completeness
          pixelval  <- matrix(initial_completeness[draw == d, initial_complete], ncol = 1)
          
          # Weighted by number of VR deaths per area in that year
          weightval <- matrix(area_child_deaths[year == y, deaths], ncol = 1)
          
          r <- 
            LogitFindK_mort(gbdval     = gbdval,
                            pixelval   = pixelval,
                            weightval  = weightval,
                            approx_0_1 = T, 
                            MaxIter    = 100)
          
          # Return as data table
          return(data.table(year = y,
                            rf = r,
                            draw = d))
        }))
      }))
    
    # Apply correction to get final completeness 
    final_child_completeness <- 
      initial_completeness %>% 
      left_join(child_logit_rf, by = "draw") %>% 
      mutate(final_complete = invlogit(logit(initial_complete) + rf)) %>% 
      data.table()
    
    # Save final completeness -- this contains the adult rf 
    saveRDS(final_child_completeness, file = paste0("<<<< FILEPATH REDACTED >>>>"))
    
    # Fit distribution to logitnormal priors pased in as arguments
    library(fitdistrplus)
    child_prior_dist <-
      rbindlist(mclapply(sort(unique(final_child_completeness$area)), mc.cores = 10, function(a) {
        message("working on area ", a)
        rbindlist(lapply(sort(unique(final_child_completeness$year)), function(y) {
          logit_norm_fit <-
            final_child_completeness %>%
            filter(area == a, year == y) %>%
            mutate(final_complete = logit(final_complete)) %>%
            pull(final_complete) %>%
            fitdist(dist = "norm", method = "mle")
          
          fit_table <- data.table(area = a,
                                  year = y,
                                  mean = logit_norm_fit$estimate[["mean"]],
                                  sd = logit_norm_fit$estimate[["sd"]])
        }))
      }))
    
    saveRDS(child_prior_dist, file = paste0("<<<< FILEPATH REDACTED >>>>"))
    
    # Make visualization plots, saved in figures folder 
    dir.create(paste0(settings_folder, "/figures/"))
    plot_completeness_prior(final_child_completeness, child_prior_dist, paste0(settings_folder, "/figures/child_prior_dist.pdf"))
    plot_dist_comparison(final_child_completeness, child_prior_dist, paste0(settings_folder, "/figures/child_prior_compare.pdf"))
  
  } else {
     child_prior_dist <- data.table(expand.grid(area = unique(initial_completeness$area), year = years))[, mean := 10][, sd := 0.01]
  }
  
  # order here is important 
  completeness_prior_dist <- 
    rbind(adult_prior_dist[, age_group := 1], child_prior_dist[, age_group := 0])[order(age_group, area, year)]
  
  saveRDS(completeness_prior_dist, file = paste0(settings_folder, "/comp_prior_dist.RDS"))
  
} else {
  message("Using completeness distribution from last model run (should not be used for model production)")
  
  # Grab saved child completeness
  if (calc_child_comp) {
    child_prior_dist <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
  } else {
    child_prior_dist <- data.table(expand.grid(area = sort(unique(shape$uid)), year = years))[, mean := 10][, sd := 0.01]
  }
  
  # Grab saved adult completeness
  if (calc_adult_comp) {
    adult_prior_dist <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
  } else {
    adult_prior_dist <- data.table(expand.grid(area = sort(unique(shape$uid)), year = years))[, mean := 10][, sd := 0.01]
  }
  
  completeness_prior_dist <- rbind(adult_prior_dist[, age_group := 1], child_prior_dist[, age_group := 0])[order(age_group, area, year)]
  
  saveRDS(completeness_prior_dist, file = paste0(settings_folder, "/comp_prior_dist.RDS"))
}

message("Done with completeness prior, now run data prep script")
######################################################
## End 
######################################################