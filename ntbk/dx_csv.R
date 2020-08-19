# clear environment and load libraries
rm(list = ls())
libs <- c('data.table', 'dplyr', 'raster', 'sf', 'fasterize')
lapply(libs, library, character.only = TRUE)

# User defined parameters
indicator_group <- 'wash'
indicator <- 'w_piped'
measure <- ''
run_date <- '2019_05_20_00_00_02/'
reg <- 'dia_south_asia'
holdout <- 0
modeling_shapefile_version <- '2019_02_27'
admin_draws <- paste0(indicator, measure,
  '_unraked_admin_draws_eb_bin0_', reg,
  '_', holdout, '.RData')
code <- "<<<< FILEPATH REDACTED >>>>"
## variable in input data that denotes ISO3 Codes
id_country_var <- 'country'
##
outdir <- "<<<< FILEPATH REDACTED >>>>"

# load data
## Loading admin_draws
outputdir <- "<<<< FILEPATH REDACTED >>>>"
setwd("<<<< FILEPATH REDACTED >>>>")

load(admin_draws)
a0_mod <- as.data.frame(admin_0)
a1_mod <- as.data.frame(admin_1)
sp_tbl <- as.data.frame(sp_hierarchy_list)

## Loading input data
setwd("<<<< FILEPATH REDACTED >>>>")
id <- list.files(pattern = paste0(indicator, '.csv'))
input_data <- as.data.frame(fread(id))
input_data$country <- input_data[,id_country_var]
input_data$prop <- input_data[,indicator]/input_data[,'N']

## Reading in shapefiles & population raster
file_dir <- "<<<< FILEPATH REDACTED >>>>"
a0_shp <- readRDS(paste0(file_dir,
                  'lbd_standard_admin_0.rds'))
a1_shp <- readRDS(paste0(file_dir,
                  'lbd_standard_admin_1.rds'))
pop_ras <- raster(
  "<<<< FILEPATH REDACTED >>>>")

# input data and admin_draws formatting function
## you can include any indicator specific data formatting in this function
## this would include any steps you perform in parallel_model that come after
## loading in input data but before running any models
format_id <- function(df, indi = indicator, repo = code, region = reg,
  shp_version = modeling_shapefile_version) {
  df <- as.data.frame(df)
  reg_list <- readRDS(paste0(repo, 'wash/00_reg_list.rds'))
  remove_col <- grep('V', names(df))
  df <- df[,setdiff(names(df), names(df)[remove_col])]
  df <- df %>%
      group_by(nid) %>%
      mutate(year = weighted.mean(x = year, w = sum_of_sample_weights)) %>%
      ungroup()
  df$N <- round(df$N, digits = 7)
  df[,indi] <- round(df[,indi], digits = 7)
  df <- filter(df, N > 0)
  df <- filter(df, year >= 2000)
  df <- filter(df, country %in% reg_list[[reg]])
  df <- filter(df, !is.na(latitude), !is.na(longitude))
  df$year <- round(df$year)

  #   # set 0s to 0.001 and 1s to 0.999
  df[,indi] <- ifelse(df[,indi]/df$N == 0, df$N*0.001, df[,indi])
  df[,indi] <- ifelse(df[,indi]/df$N == 1, df$N*0.999, df[,indi])

  # ## Remove outliers
  if (indi %in% c('s_piped', 's_imp_cr', 's_unimp_cr')) {
    vetted <- readRDS(paste0(repo, '/wash/00_sani_vetted.rds'))
    df <- filter(df, nid %in% vetted)
  } else {
    vetted <- readRDS(paste0(repo, '/wash/00_water_vetted.rds'))
    df <- filter(df, nid %in% vetted)
  }

  temp_out <- read.csv(paste0(repo, '/wash/00_temp_outlier.csv'),
  stringsAsFactors = FALSE)
  for (i in 1:nrow(temp_out)) {
    df <- filter(df, !(country == temp_out[i, 'country'] &
      year == temp_out[i, 'year']))
  }

  df <- as.data.table(df)
  return(df)
}

format_mod <- function(x) {
  geo_vars <-names(x)[grep('ADM', names(x))]
  geo_tbl <- as.data.frame(x[,geo_vars])
  names(geo_tbl) <- geo_vars

  mod_vars <-names(x)[grep('V', names(x))]
  mod_tbl <- as.data.frame(x[,mod_vars])

  geo_tbl$mean <- apply(mod_tbl, 1, mean, na.rm = TRUE)
  geo_tbl$lci <- apply(mod_tbl, 1, quantile, 0.025, na.rm = TRUE)
  geo_tbl$uci <- apply(mod_tbl, 1, quantile, 0.975, na.rm = TRUE)
  geo_tbl$year <- x$year
  return(geo_tbl)
}

# Format shps as rasters
a0_shp$ADM0_CODE <- as.numeric(as.character(a0_shp$ADM0_CODE))
a1_shp$ADM1_CODE <- as.numeric(as.character(a1_shp$ADM1_CODE))
a0_raster <- fasterize(st_as_sf(a0_shp), pop_ras, field = 'ADM0_CODE')
a1_raster <- fasterize(st_as_sf(a1_shp), pop_ras, field = 'ADM1_CODE')

# Format input data & extract adminIDs
head(format_id(input_data))
input_data <- format_id(input_data);
input_data$ADM1_CODE <- raster::extract(a1_raster, dplyr::select(input_data,
  longitude, latitude))
input_data$ADM0_CODE <- raster::extract(a0_raster, dplyr::select(input_data,
  longitude, latitude))

# Summarize input data by adminIDs
a0_input <- input_data %>%
  group_by(year, ADM0_CODE, country) %>%
  summarize(input_mean = weighted.mean(x = prop,
    w = sum_of_sample_weights*weight),
  input_ss = sum(weight*N)) %>%
  filter(!is.na(ADM0_CODE))

a1_input <- input_data %>%
  group_by(year, ADM1_CODE, ADM0_CODE, country) %>%
  summarize(input_mean = weighted.mean(x = prop,
    w = sum_of_sample_weights*weight),
  input_ss = sum(weight*N)) %>%
  filter(!is.na(ADM0_CODE))

# Format inla draws
a0_mod <- format_mod(a0_mod)
a1_mod <- format_mod(a1_mod)

# Join input data and model outputs
a0_results <- left_join(a0_mod,
  dplyr::select(a0_input, year, ADM0_CODE,
    input_ss, input_mean, country))


a1_results <- left_join(a1_mod,
  dplyr::select(a1_input, year, ADM0_CODE, ADM1_CODE,
    input_ss, input_mean))

setwd(outdir)
write.csv(a0_results, 'a0_results_dx.csv')
write.csv(a1_results, 'a0_results_dx.csv')
