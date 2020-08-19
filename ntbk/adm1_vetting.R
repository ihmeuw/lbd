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
  df <- filter(df, !is.na(latitude), !is.na(longitude))
  df$year <- round(df$year)

  # recalculate prop
  df$prop <- as.numeric(df[[indicator]])/as.numeric(df$N)

  df <- as.data.table(df)
  return(df)
}


modeling_shapefile_version <- '2019_02_27'
file_dir <- "<<<< FILEPATH REDACTED >>>>"
a0_shp <- readRDS(paste0(file_dir,
                  'lbd_standard_admin_0.rds'))
a1_shp <- readRDS(paste0(file_dir,
                  'lbd_standard_admin_1.rds'))
pop_ras <- raster("<<<< FILEPATH REDACTED >>>>")

# Format shps as rasters
a0_shp$ADM0_CODE <- as.numeric(as.character(a0_shp$ADM0_CODE))
a1_shp$ADM1_CODE <- as.numeric(as.character(a1_shp$ADM1_CODE))
a0_raster <- fasterize(st_as_sf(a0_shp), pop_ras, field = 'ADM0_CODE')
a1_raster <- fasterize(st_as_sf(a1_shp), pop_ras, field = 'ADM1_CODE')

# load data
indicator <- 'w_unimp_cr'
id <- "<<<< FILEPATH REDACTED >>>>"
input_data <- fread(id)

# Format input data
input_data <- format_id(input_data)
input_data$ADM1_CODE <- raster::extract(a1_raster, dplyr::select(input_data, longitude, latitude))

a1_input <- input_data %>%
              group_by(nid, year, ADM1_CODE, country) %>%
              summarize(mean = weighted.mean(x = prop, w = sum_of_sample_weights*weight),
                        ss = sum(weight*N))

plot_dat <- a1_input %>%
              group_by(nid, year, country) %>%
              summarize(ss = sum(ss, na.rm = TRUE), iqr = IQR(mean, na.rm = TRUE))

setwd("<<<< FILEPATH REDACTED >>>>")

pdf(paste0(indicator, '_adm1_input.pdf'),
       width = 11, height = 8.5)

  for (ccc in unique(plot_dat$country)) {
    test <- filter(plot_dat, country == ccc)
    print(
    ggplot(test) +
      geom_point(aes(x = year, y = iqr, size = ss)) +
      geom_text_repel(aes(x = year, y = iqr, label = nid)) +
      theme_bw() +
      xlim(2000, NA) +
      ylim(0, 1) +
      ggtitle(ccc)
    )
  }

dev.off()
