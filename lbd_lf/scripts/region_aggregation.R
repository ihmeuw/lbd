### Region prevalence aggregates
#############################################################

rm(list = ls())

## Set repo locations
## Get current user name
user <- Sys.info()[["user"]] ## Get current user name

## Set repo location and indicator group
core_repo <- paste0(<<<< FILEPATH REDACTED >>>>)
indic_repo <- paste0(<<<< FILEPATH REDACTED >>>>)

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- sprintf(<<<< FILEPATH REDACTED >>>>)
package_list <- c(t(read.csv(sprintf(<<<< FILEPATH REDACTED >>>>), header = FALSE)))

message("Loading in required R packages and MBG functions")
source(paste0(<<<< FILEPATH REDACTED >>>>))

mbg_setup(package_list = package_list, repos = core_repo)

## Focal 3 specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(paste0(<<<< FILEPATH REDACTED >>>>), recursive = TRUE)) {
  message(funk)
  source(paste0(<<<< FILEPATH REDACTED >>>>))
}

path <- paste0(<<<< FILEPATH REDACTED >>>>)
library(fasterize)
library(matrixStats)

######
who_pop <- fread(<<<< FILEPATH REDACTED >>>>)
who_pop <- cbind(who_pop[, 1:5], data.table(sapply(who_pop[, 6:76], function(x) {gsub(" ", "", x)})))

ig <- "lf"
ind <- "had_lf_w_resamp"
sharedir <- sprintf(<<<< FILEPATH REDACTED >>>>)

# ICT
rd_s_asia <- <<<< FILEPATH REDACTED >>>>
rd_se_asia <- <<<< FILEPATH REDACTED >>>>
rd_africa <- <<<< FILEPATH REDACTED >>>>
rd_hispaniola <- <<<< FILEPATH REDACTED >>>>
input_folder_non_mbg <- <<<< FILEPATH REDACTED >>>>

output_folder <- <<<< FILEPATH REDACTED >>>>

years <- c(1990, 1995, 2000, 2005, 2010, 2015, 2018)
year_list <- min(years):max(years)
interval_mo <- 12

thresholds <- list("_prev_0_pop_inf", "_prev_0_pop_750k", "_prev_0_pop_250k", "_prev_0_pop_100k", "_prev_0_pop_50k", "_prev_0_pop_25k", "_prev_0_pop_10k", "_prev_0_pop_5k")

method_agg <- "popweighted"

### Define regions ##########################################
searo <- c("Bangladesh", "India", "Indonesia", "Myanmar", "Nepal", "Sri Lanka", "Thailand", "Timor-Leste", "Maldives")
wpro <- c(
  "Brunei", "Cambodia", "Laos", "Malaysia", "Vietnam", "Philippines", "Papua New Guinea", "Cook Islands", "Niue", "Tonga", "Tuvalu", "Vanuatu", "Samoa", "Micronesia (Fed. States of)", "Marshall Islands", "Kiribati", "Fiji", "Palau",
  "American Samoa", "Wallis and Futuna Islands", "French Polynesia", "New Caledonia"
)
emro <- c("Egypt", "Sudan", "Yemen", "Djibouti")
amro <- c("Brazil", "Guyana", "Haiti", "Dominican Republic")
afro <- c("Angola", "Benin", "Burkina Faso", "Central African Republic", "Côte d'Ivoire", "Cameroon", "Democratic Republic of the Congo", "Republic of Congo",
          "Comoros", "Eritrea", "Ethiopia", "Gabon", "Ghana", "Guinea", "Gambia", "Guinea-Bissau", "Equatorial Guinea", "Kenya", "Liberia", "Madagascar",
          "Mali", "Mozambique", "Mauritania", "Malawi", "Niger", "Nigeria", "Senegal", "Sierra Leone", "South Sudan", "São Tomé and Príncipe", "Chad", "Togo",
          "Tanzania", "Uganda", "Zambia", "Zimbabwe")

source(<<<< FILEPATH REDACTED >>>>)
source(<<<< FILEPATH REDACTED >>>>)
source(<<<< FILEPATH REDACTED >>>>)
source(<<<< FILEPATH REDACTED >>>>)

locs <- get_location_metadata(gbd_round_id = 7, location_set_id = 35, decomp_step = "iterative")
loc_ids <- locs$location_id
gbd <- get_population(decomp_step = "iterative", gbd_round_id = 7, year_id = years, location_id = loc_ids)
gbd <- merge(as.data.table(gbd), as.data.table(locs)[, .(location_id, ihme_loc_id)])

for (p in 1:length(thresholds)) {
  pathadd <- thresholds[[p]]
  print(pathadd)
  
  ################################################################################################################################################################
  ########## Correct approach to aggregation: draws are independent between regions, so combine individual draws across regions, then compute summary stats
  ## non-MBG
  nonmbg_iso3 <- c("ASM", "BRA", "COK", "FSM", "FJI", "PYF", "GUY", "KIR", "MDV", "MHL", "NCL", "NIU", "PLW", "WSM", "TON", "TUV", "VUT", "WLF")
  geo_nonmbg <- c(
    "American Samoa", "Brazil", "Cook Islands", "Micronesia (Fed. States of)",
    "Fiji", "French Polynesia", "Guyana", "Kiribati", "Maldives", "Marshall Islands",
    "New Caledonia", "Niue", "Palau", "Samoa", "Tonga", "Tuvalu", "Vanuatu", "Wallis and Futuna Islands"
  )
  
  all_draws <- data.table()
  for (i in 1:length(nonmbg_iso3)) {
    draws <- fread(paste0(<<<< FILEPATH REDACTED >>>>))
    draws <- cbind("iso3" = nonmbg_iso3[i], "geography" = geo_nonmbg[i], draws)
    all_draws <- rbind(all_draws, draws)
  }
  
  who_pop <- as.data.frame(who_pop[who_pop$geography %in% geo_nonmbg, ])
  
  nonmbg <- copy(all_draws)
  nonmbg <- nonmbg[year %in% years]
  
  for (i in 1:nrow(nonmbg)) {
    if (nonmbg[i, iso3] %in% gbd$ihme_loc_id) {
      print(paste0("Using GBD population for ", nonmbg[i, iso3]))
      nonmbg[i, "pop"] <- gbd[ihme_loc_id == nonmbg[i, iso3] & year_id == nonmbg[i, year], population]
    } else {
      print(paste0("Using WHO population for ", nonmbg[i, iso3]))
      nonmbg[i, "pop"] <- as.numeric(who_pop[who_pop$geography == nonmbg[i, geography], which(colnames(who_pop) == nonmbg[i, year])]) * 1000
    }
  }
  
  # WHO populations used for: French Polynesia, New Caledonia, Wallis and Futuna
  
  ## S/SE Asia MBG models
  input_dir_s_asia <- paste0(<<<< FILEPATH REDACTED >>>>)
  output_dir_s_asia <- paste0(<<<< FILEPATH REDACTED >>>>)
  load(paste0(<<<< FILEPATH REDACTED >>>>))
  adm0_agg_s_asia <- fread(paste0(<<<< FILEPATH REDACTED >>>>))
  admin_0_s_asia <- copy(admin_0)
  
  input_dir_se_asia <- paste0(<<<< FILEPATH REDACTED >>>>)
  output_dir_se_asia <- paste0(<<<< FILEPATH REDACTED >>>>)
  load(paste0(<<<< FILEPATH REDACTED >>>>))
  adm0_agg_se_asia <- fread(paste0(<<<< FILEPATH REDACTED >>>>))
  admin_0_se_asia <- copy(admin_0)
  
  adm0_agg_asia <- rbind(adm0_agg_s_asia, adm0_agg_se_asia)
  admin_0 <- rbind(admin_0_s_asia, admin_0_se_asia)
  
  admin_0[, c("gam", "gbm", "lasso") := NULL]
  
  admin_0 <- merge(admin_0, unique(adm0_agg_asia[, .(ADM0_CODE, ADM0_NAME)]), by = "ADM0_CODE", all.x = TRUE)
  
  nonmbg[, "iso3" := NULL]
  admin_0 <- cbind("geography" = admin_0$ADM0_NAME, admin_0[, 2:1003])
  
  all_draws <- rbind(nonmbg, copy(admin_0))
  
  #### Hispaniola draws
  input_dir_hispaniola <- paste0(<<<< FILEPATH REDACTED >>>>)
  output_dir_hispaniola <- paste0(<<<< FILEPATH REDACTED >>>>)
  load(paste0(<<<< FILEPATH REDACTED >>>>))
  adm0_agg_hispaniola <- fread(paste0(<<<< FILEPATH REDACTED >>>>))
  
  admin_0[, c("gam", "gbm", "lasso") := NULL]
  
  admin_0 <- merge(admin_0, unique(adm0_agg_hispaniola[, .(ADM0_CODE, ADM0_NAME)]), by = "ADM0_CODE", all.x = TRUE)
  admin_0 <- cbind("geography" = admin_0$ADM0_NAME, admin_0[, 2:1003])
  
  all_draws <- rbind(all_draws, copy(admin_0))
  
  #### Africa draws
  input_dir_africa <- paste0(<<<< FILEPATH REDACTED >>>>)
  output_dir_africa <- paste0(<<<< FILEPATH REDACTED >>>>)
  load(paste0(<<<< FILEPATH REDACTED >>>>))
  adm0_agg_africa <- fread(paste0(<<<< FILEPATH REDACTED >>>>))
  
  admin_0[, c("gam", "gbm", "lasso") := NULL]
  
  admin_0 <- merge(admin_0, unique(adm0_agg_africa[, .(ADM0_CODE, ADM0_NAME)]), by = "ADM0_CODE", all.x = TRUE)
  admin_0 <- cbind("geography" = admin_0$ADM0_NAME, admin_0[, 2:1003])
  
  all_draws <- rbind(all_draws, copy(admin_0))
  
  #### Subset to years
  admin_0 <- all_draws[year %in% years]

  locs[location_name  == "Brunei Darussalam", location_name := "Brunei"]
  locs[location_name  == "Lao People's Democratic Republic", location_name := "Laos"]
  locs[location_name  == "Micronesia (Federated States of)", location_name := "Micronesia (Fed. States of)"]
  locs[location_name  == "Congo", location_name := "Republic of Congo"]
  locs[location_name  == "Sao Tome and Principe", location_name := "São Tomé and Príncipe"]
  locs[location_name  == "United Republic of Tanzania", location_name := "Tanzania"]
  locs[location_name  == "Viet Nam", location_name := "Vietnam"]
  
  
  admin_0 <- merge(admin_0, locs[(location_type == "admin0") | (ihme_loc_id %in% c("ASM", "COK", "NIU")), c("location_name", "ihme_loc_id")], by.x = "geography", by.y = "location_name", all.x = TRUE)
  
  for (i in 1:length(admin_0)) {
    if (admin_0[i, ihme_loc_id] %in% gbd$ihme_loc_id) {
      admin_0[i, "pop_gbd"] <- gbd[ihme_loc_id == admin_0[i, ihme_loc_id] & year_id == admin_0[i, year], population]
    }
  }
  
  admin_0[!is.na(pop_gbd), "final_pop" := pop_gbd]
  admin_0[is.na(pop_gbd), "final_pop" := pop]
  
  #### Now perform regional aggregations
  searo_cases <- lapply(years, FUN = function(i, method = method_agg) {
    if (method == "popweighted") { # Method with national pop-weighted prevalence. Move to case space first before calculating summary stats in prev space.
      prevs <- admin_0[year == i, grep("V", names(admin_0)), with = F]
      pops <- admin_0[year == i, final_pop]
      cases <- apply(prevs, 2, FUN = function(i) i * pops)
      lower_cases <- as.matrix(cases) %>% rowQuantiles(probs = .025)
      upper_cases <- as.matrix(cases) %>% rowQuantiles(probs = .975)
      mean_cases <- as.matrix(cases) %>% rowMeans()
      temp <- as.data.table(cbind("geography" = admin_0[year == i, geography], lower_cases, mean_cases, upper_cases))
      temp <- temp[geography %in% searo, ]
      country_out <- copy(temp)
      country_out[, year := i]
      
      # region aggregation
      subset <- admin_0[year == i & geography %in% searo, ]
      prevs <- subset[year == i, grep("V", names(admin_0)), with = F]
      pops <- subset[year == i, final_pop]
      cases <- apply(prevs, 2, FUN = function(i) sum(i * pops))
      lower_cases <- quantile(cases, probs = .025)
      upper_cases <- quantile(cases, probs = .975)
      mean_cases <- mean(cases)
      temp <- c(lower_cases, mean_cases, upper_cases)
    }
    
    return(list(temp, country_out))
  })
  
  country_cases <- rbindlist(lapply(searo_cases, FUN = function(i) i[[2]]), fill = T, use.names = T)
  searo_cases <- lapply(searo_cases, FUN = function(i) i[[1]])
  searo_cases <- do.call(rbind, searo_cases)
  colnames(searo_cases)[2] <- "mean_cases"
  
  
  wpro_cases <- lapply(years, FUN = function(i, method = method_agg) {
    if (method == "popweighted") { # Method with national pop-weighted prevalence. Move to case space first before calculating summary stats in prev space.
      prevs <- admin_0[year == i, grep("V", names(admin_0)), with = F]
      pops <- admin_0[year == i, final_pop]
      cases <- apply(prevs, 2, FUN = function(i) i * pops)
      lower_cases <- as.matrix(cases) %>% rowQuantiles(probs = .025)
      upper_cases <- as.matrix(cases) %>% rowQuantiles(probs = .975)
      mean_cases <- as.matrix(cases) %>% rowMeans()
      temp <- as.data.table(cbind("geography" = admin_0[year == i, geography], lower_cases, mean_cases, upper_cases))
      temp <- temp[geography %in% wpro, ]
      country_out <- copy(temp)
      country_out[, year := i]
      
      # region aggregation
      subset <- admin_0[year == i & geography %in% wpro, ]
      prevs <- subset[year == i, grep("V", names(admin_0)), with = F]
      pops <- subset[year == i, final_pop]
      cases <- apply(prevs, 2, FUN = function(i) sum(i * pops))
      lower_cases <- quantile(cases, probs = .025)
      upper_cases <- quantile(cases, probs = .975)
      mean_cases <- mean(cases)
      temp <- c(lower_cases, mean_cases, upper_cases)
    }
    
    return(list(temp, country_out))
  })
  
  country_cases <- rbind(country_cases, rbindlist(lapply(wpro_cases, FUN = function(i) i[[2]]), fill = T, use.names = T))
  wpro_cases <- lapply(wpro_cases, FUN = function(i) i[[1]])
  wpro_cases <- do.call(rbind, wpro_cases)
  colnames(wpro_cases)[2] <- "mean_cases"
  
  
  amro_cases <- lapply(years, FUN = function(i, method = method_agg) {
    if (method == "popweighted") { # Method with national pop-weighted prevalence. Move to case space first before calculating summary stats in prev space.
      prevs <- admin_0[year == i, grep("V", names(admin_0)), with = F]
      pops <- admin_0[year == i, final_pop]
      cases <- apply(prevs, 2, FUN = function(i) i * pops)
      lower_cases <- as.matrix(cases) %>% rowQuantiles(probs = .025)
      upper_cases <- as.matrix(cases) %>% rowQuantiles(probs = .975)
      mean_cases <- as.matrix(cases) %>% rowMeans()
      temp <- as.data.table(cbind("geography" = admin_0[year == i, geography], lower_cases, mean_cases, upper_cases))
      temp <- temp[geography %in% amro, ]
      country_out <- copy(temp)
      country_out[, year := i]
      
      # region aggregation
      subset <- admin_0[year == i & geography %in% amro, ]
      prevs <- subset[year == i, grep("V", names(admin_0)), with = F]
      pops <- subset[year == i, final_pop]
      cases <- apply(prevs, 2, FUN = function(i) sum(i * pops))
      lower_cases <- quantile(cases, probs = .025)
      upper_cases <- quantile(cases, probs = .975)
      mean_cases <- mean(cases)
      temp <- c(lower_cases, mean_cases, upper_cases)
    }
    
    return(list(temp, country_out))
  })
  
  country_cases <- rbind(country_cases, rbindlist(lapply(amro_cases, FUN = function(i) i[[2]]), fill = T, use.names = T))
  amro_cases <- lapply(amro_cases, FUN = function(i) i[[1]])
  amro_cases <- do.call(rbind, amro_cases)
  colnames(amro_cases)[2] <- "mean_cases"
  
  
  emro_cases <- lapply(years, FUN = function(i, method = method_agg) {
    if (method == "popweighted") { # Method with national pop-weighted prevalence. Move to case space first before calculating summary stats in prev space.
      prevs <- admin_0[year == i, grep("V", names(admin_0)), with = F]
      pops <- admin_0[year == i, final_pop]
      cases <- apply(prevs, 2, FUN = function(i) i * pops)
      lower_cases <- as.matrix(cases) %>% rowQuantiles(probs = .025)
      upper_cases <- as.matrix(cases) %>% rowQuantiles(probs = .975)
      mean_cases <- as.matrix(cases) %>% rowMeans()
      temp <- as.data.table(cbind("geography" = admin_0[year == i, geography], lower_cases, mean_cases, upper_cases))
      temp <- temp[geography %in% emro, ]
      country_out <- copy(temp)
      country_out[, year := i]
      
      # region aggregation
      subset <- admin_0[year == i & geography %in% emro, ]
      prevs <- subset[year == i, grep("V", names(admin_0)), with = F]
      pops <- subset[year == i, final_pop]
      cases <- apply(prevs, 2, FUN = function(i) sum(i * pops))
      lower_cases <- quantile(cases, probs = .025)
      upper_cases <- quantile(cases, probs = .975)
      mean_cases <- mean(cases)
      temp <- c(lower_cases, mean_cases, upper_cases)
    }
    
    return(list(temp, country_out))
  })
  
  country_cases <- rbind(country_cases, rbindlist(lapply(emro_cases, FUN = function(i) i[[2]]), fill = T, use.names = T))
  emro_cases <- lapply(emro_cases, FUN = function(i) i[[1]])
  emro_cases <- do.call(rbind, emro_cases)
  colnames(emro_cases)[2] <- "mean_cases"
  
  
  afro_cases <- lapply(years, FUN = function(i, method = method_agg) {
    if (method == "popweighted") { # Method with national pop-weighted prevalence. Move to case space first before calculating summary stats in prev space.
      prevs <- admin_0[year == i, grep("V", names(admin_0)), with = F]
      pops <- admin_0[year == i, final_pop]
      cases <- apply(prevs, 2, FUN = function(i) i * pops)
      lower_cases <- as.matrix(cases) %>% rowQuantiles(probs = .025)
      upper_cases <- as.matrix(cases) %>% rowQuantiles(probs = .975)
      mean_cases <- as.matrix(cases) %>% rowMeans()
      temp <- as.data.table(cbind("geography" = admin_0[year == i, geography], lower_cases, mean_cases, upper_cases))
      temp <- temp[geography %in% afro, ]
      country_out <- copy(temp)
      country_out[, year := i]
      
      # region aggregation
      subset <- admin_0[year == i & geography %in% afro, ]
      prevs <- subset[year == i, grep("V", names(admin_0)), with = F]
      pops <- subset[year == i, final_pop]
      cases <- apply(prevs, 2, FUN = function(i) sum(i * pops))
      lower_cases <- quantile(cases, probs = .025)
      upper_cases <- quantile(cases, probs = .975)
      mean_cases <- mean(cases)
      temp <- c(lower_cases, mean_cases, upper_cases)
    }
    
    return(list(temp, country_out))
  })
  
  country_cases <- rbind(country_cases, rbindlist(lapply(afro_cases, FUN = function(i) i[[2]]), fill = T, use.names = T))
  afro_cases <- lapply(afro_cases, FUN = function(i) i[[1]])
  afro_cases <- do.call(rbind, afro_cases)
  colnames(afro_cases)[2] <- "mean_cases"
  
  
  ### Regional totals #######################################################################################
  all_regions <- list(afro_cases, amro_cases, emro_cases, searo_cases, wpro_cases)
  names(all_regions) <- c("afro", "amro", "emro", "searo", "wpro")
  for (i in 1:length(all_regions)) {
    all_regions[[i]] <- cbind(as.data.table(all_regions[[i]]), years)
    all_regions[[i]][, reg := names(all_regions)[i]]
  }
  all_regions <- rbindlist(all_regions)
  
  all_regions[, 1:3] <- round(all_regions[, 1:3], 0)
  country_cases$lower_cases <- round(as.numeric(country_cases$lower_cases), 0)
  country_cases$mean_cases <- round(as.numeric(country_cases$mean_cases), 0)
  country_cases$upper_cases <- round(as.numeric(country_cases$upper_cases), 0)
  
  write.csv(all_regions, paste0(<<<< FILEPATH REDACTED >>>>), row.names = F)
  write.csv(country_cases, paste0(<<<< FILEPATH REDACTED >>>>), row.names = F)
  
  
  ### Global totals #######################################################################################
  global_cases <- lapply(years, FUN = function(i, method = method_agg) {
    if (method == "popweighted") { # Method with national pop-weighted prevalence. Move to case space first before calculating summary stats in prev space.
      subset <- admin_0[year == i, ]
      prevs <- subset[year == i, grep("V", names(admin_0)), with = F]
      pops <- subset[year == i, final_pop]
      cases <- apply(prevs, 2, FUN = function(i) sum(i * pops))
      lower_cases <- quantile(cases, probs = .025)
      upper_cases <- quantile(cases, probs = .975)
      mean_cases <- mean(cases)
      temp <- c(lower_cases, mean_cases, upper_cases)
    }
    
    return(list(temp))
  })
  
  global_cases <- lapply(global_cases, FUN = function(i) i[[1]])
  global_cases <- do.call(rbind, global_cases)
  colnames(global_cases)[2] <- "mean_cases"
  
  write.csv(cbind("year"=years, round(global_cases, 0)), paste0(<<<< FILEPATH REDACTED >>>>), row.names = F)
  write.csv(admin_0[, c("geography", "year", "ihme_loc_id", "final_pop")], paste0(<<<< FILEPATH REDACTED >>>>), row.names = F)
}
