# Clear environment
rm(list = ls())

# Set up environment with packages and lbd_core functions
user <- Sys.info()['user']
repo <- "<<<< FILEPATH REDACTED >>>>"
commondir <- "<<<< FILEPATH REDACTED >>>>"
package_list <- c(t(read.csv(sprintf("%s/package_list.csv", commondir), header = FALSE)))
setwd(repo)
core_repo <- repo
libs <- c('seegSDM', 'seegMBG', 'data.table',
  'mgcv', 'raster', 'sf', 'fasterize', 'dplyr', 'ggplot2', package_list)
lapply(libs, library, character.only = TRUE)
message("Loading in required R packages and MBG functions")
source(paste0(repo, "/mbg_central/setup.R"))
mbg_setup(package_list = package_list, repos = repo)

# User defined variables
indi <- 'w_piped'
reg <- 'dia_wssa-civ'
run_date <- '2019_09_29_00_00_01'
imprvd <- 'TRUE'
# User defined functions
format_id <- function(df, indi = indicator) {
  df <- as.data.frame(df)
  remove_col <- grep('V', names(df))
  df <- df[,setdiff(names(df), names(df)[remove_col])]
  user            <- Sys.info()['user']
  repo            <- "<<<< FILEPATH REDACTED >>>>"
  reg_list <- readRDS(paste0(repo, "wash/00_reg_list.rds"))

  df <- dplyr::filter(df, country %in% reg_list[['dia_wssa-civ']])
  df <- df %>%
      group_by(nid) %>%
      mutate(year = weighted.mean(x = year, w = sum_of_sample_weights, na.rm = TRUE)) %>%
      ungroup()
  df$N <- round(df$N, digits = 7)
  df[,indi] <- round(df[,indi], digits = 7)
  df <- dplyr::filter(df, N > 0)
  df <- dplyr::filter(df, year >= 2000)
  df <- dplyr::filter(df, !is.na(latitude), !is.na(longitude))
  df$year <- round(df$year)

  print(nrow(df))
  # recalculate prop
	repo <- "<<<< FILEPATH REDACTED >>>>"
  ## Remove outliers
  if (indi %in% c('s_piped', 's_imp_cr', 's_unimp_cr', 's_od')) {
    vetted <- readRDS(paste0(repo, '/wash/00_sani_vetted.rds'))
    df <- dplyr::filter(df, nid %in% vetted)
  } else {
    vetted <- readRDS(paste0(repo, '/wash/00_water_vetted.rds'))
    df <- dplyr::filter(df, nid %in% vetted)
  }

  print(nrow(df))
  df$survey_series[grep('JHSPH', df$survey_series)] <- 'JHSPH'
  return(df)
}
no_shp_factors <- function(df) {
  codes <- grep('CODE', names(df))
  for (i in codes) {
    df[,i] <- as.numeric(as.character(df[,i]))
  }

  names <- grep('NAME', names(df))
  for (i in names) {
    df[,i] <- as.character(df[,i])
  }
  return(df)
}
# Load data
## Spatial objects
### Load Shapefiles
file_dir <- "<<<< FILEPATH REDACTED >>>>"
a0_shp <- readRDS(paste0(file_dir,
                  'lbd_standard_admin_0.rds'))
a0_shp@data <- no_shp_factors(a0_shp@data)
a1_shp <- readRDS(paste0(file_dir,
                  'lbd_standard_admin_1.rds'))
a1_shp@data <- no_shp_factors(a1_shp@data)
a2_shp <- readRDS(paste0(file_dir,
                  'lbd_standard_admin_2.rds'))
a2_shp@data <- no_shp_factors(a2_shp@data)

### Load pop raster
pop_ras <- raster("<<<< FILEPATH REDACTED >>>>")
## Model draws
for (iii in 0:2) {
  setwd("<<<< FILEPATH REDACTED >>>>")
  load(paste0(indi, '_unraked_admin_draws_eb_bin0_', reg,'_0.RData'))
  piped <- as.data.frame(get(paste0('admin_', iii)))
  pdraw <- as.matrix(as.data.frame(piped)[,grep('V', names(piped))])

  if (imprvd) {
    setwd("<<<< FILEPATH REDACTED >>>>")
    load(paste0('w_imp_unraked_admin_draws_eb_bin0_', reg,'_0.RData'))
    imp_cr <- as.data.frame(get(paste0('admin_', iii)))
    idraw <- as.matrix(as.data.frame(imp_cr)[,grep('V', names(imp_cr))])
    impdraw <- idraw
  } else {
    impdraw <- pdraw
  }

  assign(paste0('piped_',iii), piped)
  assign(paste0('impdraw_',iii), impdraw)
}

## Input data for model
pdat <- read.csv("<<<< FILEPATH REDACTED >>>>", stringsAsFactors = FALSE)
idat <- read.csv("<<<< FILEPATH REDACTED >>>>", stringsAsFactors = FALSE)

# Process data
## Make admin rasters
a0_shp <- a0_shp[which(a0_shp@data$ADM0_NAME == 'Nigeria'),]
a1_shp <- a1_shp[which(a1_shp@data$ADM0_NAME == 'Nigeria'),]
a2_shp <- a2_shp[which(a2_shp@data$ADM0_NAME == 'Nigeria'),]

a0_raster <- fasterize(st_as_sf(a0_shp), pop_ras, field = 'ADM0_CODE')
a1_raster <- fasterize(st_as_sf(a1_shp), pop_ras, field = 'ADM1_CODE')
a2_raster <- fasterize(st_as_sf(a2_shp), pop_ras, field = 'ADM2_CODE')

## Model output
for (iii in 0:2) {
  impdraw <- get(paste0('impdraw_',iii))
  improved <- as.data.frame(cbind(
    apply(impdraw, 1, mean, na.rm = TRUE),
    apply(impdraw, 1, quantile, probs = 0.025,  na.rm = TRUE),
    apply(impdraw, 1, quantile, probs = 0.975, na.rm = TRUE)))
  names(improved)[1] <- 'mean'
  names(improved)[2] <- 'lower'
  names(improved)[3] <- 'upper'
  print(head(cbind(get(paste0('piped_',iii))[,1:2], improved)))
  improved <- cbind(get(paste0('piped_',iii))[,1:2], improved)
  ind_imp <- improved[which(improved[,paste0('ADM', iii, '_CODE')] %in% get(paste0('a', iii, '_shp'))@data[,paste0('ADM', iii, '_CODE')]),]
  ind_imp <- left_join(ind_imp, distinct(dplyr::select(get(paste0('a', iii, '_shp'))@data, paste0('ADM', iii, '_CODE'),
    paste0('ADM', iii, '_NAME'))))
  assign(paste0('ind_imp_',iii), ind_imp)
}

## Model input
if (indi == 's_piped') {
  pdat$year <- pdat$surv_year
}
pdat <- format_id(pdat, indi = indi)
pdat$ADM0_CODE <- raster::extract(a0_raster,
  dplyr::select(pdat, longitude, latitude))
pdat$ADM1_CODE <- raster::extract(a1_raster,
  dplyr::select(pdat, longitude, latitude))
pdat$ADM2_CODE <- raster::extract(a2_raster,
  dplyr::select(pdat, longitude, latitude))
pdat <- as.data.frame(pdat)
pdat_0 <- dplyr::filter(pdat, country == 'NGA') %>%
  dplyr::group_by(nid, year_median, ADM0_CODE, survey_series, point) %>%
  dplyr::summarize(prev = weighted.mean(x = prop, w = weight*N), ss = sum(weight*N))
pdat_1 <- dplyr::filter(pdat, country == 'NGA') %>%
  dplyr::group_by(nid, year_median, ADM1_CODE, survey_series, point) %>%
  dplyr::summarize(prev = weighted.mean(x = prop, w = weight*N), ss = sum(weight*N))
pdat_2 <- dplyr::filter(pdat, country == 'NGA') %>%
  dplyr::group_by(nid, year_median, ADM2_CODE, survey_series, point) %>%
  dplyr::summarize(prev = weighted.mean(x = prop, w = weight*N), ss = sum(weight*N))

if (imprvd) {
  idat <- format_id(idat, indi = 'w_imp_cr')
  idat$ADM0_CODE <- raster::extract(a0_raster,
    dplyr::select(idat, longitude, latitude))
  idat$ADM1_CODE <- raster::extract(a1_raster,
    dplyr::select(idat, longitude, latitude))
  idat$ADM2_CODE <- raster::extract(a2_raster,
    dplyr::select(idat, longitude, latitude))

    idat_0 <- dplyr::filter(idat, country == 'NGA') %>%
      dplyr::group_by(nid, year_median, ADM0_CODE, survey_series, point) %>%
      dplyr::summarize(prev = weighted.mean(x = prop, w = weight*N))
    idat_1 <- dplyr::filter(idat, country == 'NGA') %>%
      dplyr::group_by(nid, year_median, ADM1_CODE, survey_series, point) %>%
      dplyr::summarize(prev = weighted.mean(x = prop, w = weight*N))
    idat_2 <- dplyr::filter(idat, country == 'NGA') %>%
      dplyr::group_by(nid, year_median, ADM2_CODE, survey_series, point) %>%
      dplyr::summarize(prev = weighted.mean(x = prop, w = weight*N))

  mdat_0 <- left_join(pdat_0, dplyr::rename(idat_0, imp_cr = prev))
  mdat_1 <- left_join(pdat_1, dplyr::rename(idat_1, imp_cr = prev))
  mdat_2 <- left_join(pdat_2, dplyr::rename(idat_2, imp_cr = prev))

  mdat_0 <- mdat_0[complete.cases(mdat_0),] %>%
    mutate(prev = prev+((1 - prev)*imp_cr)) %>%
    dplyr::select(-imp_cr)
  mdat_1 <- mdat_1[complete.cases(mdat_1),] %>%
    mutate(prev = prev+((1 - prev)*imp_cr)) %>%
    dplyr::select(-imp_cr)
  mdat_2 <- mdat_2[complete.cases(mdat_2),] %>%
    mutate(prev = prev+((1 - prev)*imp_cr)) %>%
    dplyr::select(-imp_cr)

} else {
  mdat_0 <- pdat_0
  mdat_1 <- pdat_1
  mdat_2 <- pdat_2
}

mdat_0 <- left_join(mdat_0, distinct(dplyr::select(a0_shp@data, 'ADM0_CODE', 'ADM0_NAME')))
mdat_1 <- left_join(mdat_1, distinct(dplyr::select(a1_shp@data, 'ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME')))
mdat_2 <- left_join(mdat_2, distinct(dplyr::select(a2_shp@data, 'ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME',
  'ADM2_CODE', 'ADM2_NAME')))

mdat_0 <- dplyr::rename(mdat_0, year = year_median, svy_id = nid, source = survey_series, outcome = prev, N = ss) %>%
  dplyr::filter(!is.na(ADM0_CODE)) %>% dplyr::distinct()
mdat_1 <- dplyr::rename(mdat_1, year = year_median, svy_id = nid, source = survey_series, outcome = prev, N = ss)  %>%
  dplyr::filter(!is.na(ADM1_CODE)) %>% dplyr::distinct()
mdat_2 <- dplyr::rename(mdat_2, year = year_median, svy_id = nid, source = survey_series, outcome = prev, N = ss)  %>%
  dplyr::filter(!is.na(ADM2_CODE)) %>% dplyr::distinct()
# Save plots
## Fast Mosser
pdf("<<<< FILEPATH REDACTED >>>>", width = 11, height = 8)
print(
ggplot() +
	geom_line(data = ind_imp, aes(x = year, y = prev)) +
  geom_ribbon(data = ind_imp, aes(x = year, ymin = lci, ymax = uci), alpha = 0.3) +
	geom_point(data = mdat, aes(x = year, y = prev, col = survey_series,
		size = ss), alpha = 0.8) +
	facet_wrap(. ~ ADM1_NAME) +
	theme_bw() +
  ggtitle('Piped Water')
)
dev.off()
ind_imp_0$mean <- 1
ind_imp_1$mean <- 1
ind_imp_2$mean <- 1

## Full Mosser
subnational_ts_plots(ad0_df = ind_imp_0,
                     ad2_df = ind_imp_2,
                     ad1_df = ind_imp_1,
                     ad0_data = mdat_0,
                     ad1_data = mdat_1,
                     ad2_data = mdat_2,
                     ad0_shape = a0_shp,
                     ad1_shape = a1_shp,
                     ad2_shape = a2_shp,
                     ind_title = 'Improved Water',
                     highisbad = F,
                     ad0_map_regions = 'nga',
                     verbose = T,
                     plot_data = T,
                     plot_levels = c('ad1', 'ad2'),
                     out_dir = "<<<< FILEPATH REDACTED >>>>",
                     val_range = c(0,1))

# Write out processed data
setwd("<<<< FILEPATH REDACTED >>>>")
write.csv(ind_imp, 'water_a1_mod_agg.csv', row.names = FALSE)
write.csv(mdat, 'water_a1_idat_agg.csv', row.names = FALSE)
