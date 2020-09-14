#####################################################################
# Make IQR plot by admin/NID for input data
#####################################################################

# (1) Setup -------------------------------------------------------------------------------------
rm(list = ls())

#user inputs
indicator <- 'has_lri'
indicator_group <- 'lri'
data_tag <- '' #tag for input data set
modeling_shapefile_version <- '2019_09_10'
input_data_date <- '2020_04_08' #date collapse was run
out_dir <- '<<<< FILEPATH REDACTED >>>>'
dir.create(out_dir, recursive = TRUE)

# load libraries
libs <- c('ggplot2', 'data.table', 'dplyr', 'ggpubr', 'raster', 'sf', 'fasterize', 'ggrepel')
lapply(libs, library, character.only = TRUE)
select <- dplyr::select

#define function to format data
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

# (2) Load data and load shapefiles, pop ---------------------------------------------------------

id <- '<<<< FILEPATH REDACTED >>>>'
input_data <- fread(id)

file_dir <- '<<<< FILEPATH REDACTED >>>>'
a0_shp <- readRDS(paste0(file_dir,
                         'lbd_standard_admin_0.rds'))
a1_shp <- readRDS(paste0(file_dir,
                         'lbd_standard_admin_1.rds'))
pop_ras <- raster('<<<< FILEPATH REDACTED >>>>')

# Format shps as rasters
a0_shp$ADM0_CODE <- as.numeric(as.character(a0_shp$ADM0_CODE))
a1_shp$ADM1_CODE <- as.numeric(as.character(a1_shp$ADM1_CODE))
a0_raster <- fasterize(st_as_sf(a0_shp), pop_ras, field = 'ADM0_CODE')
a1_raster <- fasterize(st_as_sf(a1_shp), pop_ras, field = 'ADM1_CODE')

# Format input data
input_data <- format_id(input_data)
input_data$ADM1_CODE <- raster::extract(a1_raster, dplyr::select(input_data, longitude, latitude))
a1_input <- input_data %>%
  group_by(nid, year, ADM1_CODE, country, survey_series) %>%
  summarize(mean = weighted.mean(x = prop, w = sum_of_sample_weights*weight),
            ss = sum(weight*N))
plot_dat <- a1_input %>%
  group_by(nid, year, country, survey_series) %>%
  summarize(ss = sum(ss, na.rm = TRUE), iqr = IQR(mean, na.rm = TRUE), n_admins_surveyed = length(unique(ADM1_CODE)))

iso <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  select(iso3, ADM0_NAME)

plot_dat <- merge(plot_dat, iso, by.x = 'country', by.y = 'iso3')

# (3) Plot data ------------------------------------------------------------------------------

setwd(out_dir)
write.csv(a1_input, paste0(indicator, '_a1_input.csv'))
pdf(paste0(indicator, '_adm1_input.pdf'),
    width = 11, height = 8.5)

max_iqr <- max(plot_dat$iqr)
for (ccc in unique(plot_dat$country)) {
  test <- filter(plot_dat, country == ccc)
  print(
    ggplot(test) +
      geom_point(aes(x = year, y = iqr, size = ss, color = survey_series)) +
      geom_text_repel(aes(x = year, y = iqr, label = nid), box.padding = 1) +
      geom_label_repel(aes(x = year, y = iqr, label = n_admins_surveyed)) +
      theme_bw() + 
      xlim(2000, NA) +
      ylim(0, max_iqr) +
      ggtitle(paste(ccc, unique(test$ADM0_NAME)))
  )   
}
dev.off()

