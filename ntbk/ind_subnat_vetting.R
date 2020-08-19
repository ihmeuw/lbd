# load libraries
libs <- c('ggplot2', 'data.table', 'dplyr', 'ggpubr', 'raster', 'sf', 'fasterize', 'ggrepel')
lapply(libs, library, character.only = TRUE)

# format data
format_id <- function(df, indi = indicator) {
  df <- as.data.frame(df)
  remove_col <- grep('V', names(df))
  df <- df[,setdiff(names(df), names(df)[remove_col])]
  df <- df %>%
      group_by(nid) %>%
      mutate(year = weighted.mean(x = year, w = sum_of_sample_weights)) %>%
      ungroup()
  df$N <- round(df$N, digits = 7)
  df[,indicator] <- round(df[,indicator], digits = 7)
  df <- filter(df, N > 0)
  df <- filter(df, year >= 2000)
  df <- filter(df, country %in% reg_list[[reg]])
  df <- filter(df, !is.na(latitude), !is.na(longitude))
  df$year <- round(df$year)

    # set 0s to 0.001 and 1s to 0.999
  df[,indicator] <- ifelse(df[,indicator]/df$N == 0, df$N*0.001, df[,indicator])
  df[,indicator] <- ifelse(df[,indicator]/df$N == 1, df$N*0.999, df[,indicator])

  # recalculate prop
  df$prop <- as.numeric(df[[indicator]]/df$N)

  df <- as.data.table(df)
  return(df)
}

# load data
modeling_shapefile_version <- '2019_09_10'
file_dir <- "<<<< FILEPATH REDACTED >>>>"
a0_shp <- readRDS(paste0(file_dir,
                  'lbd_standard_admin_0.rds'))
a1_shp <- readRDS(paste0(file_dir,
                  'lbd_standard_admin_1.rds'))
pop_ras <- raster("<<<< FILEPATH REDACTED >>>>")
reg_list <- list('ind_pak' = c('IND', 'PAK'))

indicator <- 's_network_cr'
input_data <- fread("<<<< FILEPATH REDACTED >>>>")

# Format shps as rasters
a0_shp$ADM0_CODE <- as.numeric(as.character(a0_shp$ADM0_CODE))
a1_shp$ADM1_CODE <- as.numeric(as.character(a1_shp$ADM1_CODE))
a0_raster <- fasterize(st_as_sf(a0_shp), pop_ras, field = 'ADM0_CODE')
a1_raster <- fasterize(st_as_sf(a1_shp), pop_ras, field = 'ADM1_CODE')

# Format input data
input_data <- filter(input_data, country == 'IND')
input_data$ADM0_CODE <- 105
input_data$ADM1_CODE <- raster::extract(a1_raster, dplyr::select(input_data, longitude, latitude))

a0_input <- input_data %>%
              group_by(nid, year, ADM0_CODE) %>%
              summarize(mean = weighted.mean(x = prop, w = sum_of_sample_weights*weight),
                        ss = sum(weight*N))

a1_input <- input_data %>%
              group_by(nid, year, ADM1_CODE, ADM0_CODE) %>%
              summarize(mean = weighted.mean(x = prop, w = sum_of_sample_weights*weight),
                        ss = sum(weight*N))
a1_sp_dat <- a1_shp@data
a1_sp_dat$ADM1_CODE <- as.numeric(as.character(a1_sp_dat$ADM1_CODE))
a1_sp_dat$ADM1_NAME <- as.character(a1_sp_dat$ADM1_NAME)
a1_sp_dat <- dplyr::select(a1_sp_dat, ADM1_CODE, ADM1_NAME)
a1_input <- left_join(a1_input, a1_sp_dat)

pdf("<<<< FILEPATH REDACTED >>>>")
for (i in unique(a1_input$ADM1_NAME)) {
  print(i)
  if (!is.na(i)) {
    plotdat <- filter(a1_input, ADM1_NAME == i)
    gg1 <- ggplot() +
      geom_point(aes(x = year, y = mean, size = ss), data = plotdat) +
      geom_text_repel(aes(x = year, y = mean, label = nid), data = plotdat) +
      xlim(2000, NA) + ylim(0, 1) +
      theme_bw() +
      ggtitle(i)
    print(gg1)
  }
}
dev.off()
