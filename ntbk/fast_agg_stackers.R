
rundate <- '2019_09_01_00_00_00'
region <- 'dia_south_asia'
indicator <- 'w_piped'

setwd("<<<< FILEPATH REDACTED >>>>")

source('<<<< FILEPATH REDACTED >>>>/mbg/post_est_custom/grab_fx.R')
source('<<<< FILEPATH REDACTED >>>>/mbg/post_est_custom/raster_postest.R')

admin_lvl <- 1
stackers <- grab_stacker_raster('w_piped',
  '2019_09_01_00_00_00',
  'dia_south_asia',
  'wash',
  0)

message("Loading shapefile")
file_dir <- "<<<< FILEPATH REDACTED >>>>"
shapes <- readRDS(paste0(file_dir, 'lbd_standard_admin_1.rds'))
shapes2 <- shapes[which(shapes$ADM0_NAME == 'India'),]
file_dir <- "<<<< FILEPATH REDACTED >>>>"
pop <- raster(paste0(file_dir, 'worldpop_total_1y_', '2010', '_00_00.tif'))
pop <- crop(pop, stackers[[1]][[1]])

stack_agg <- do.call(rbind, lapply(1:length(stackers), function(y) {
  results <- do.call(rbind, lapply(1:18, function(x, z = y) {
    tbl <- agg_raster(shp = shapes2, weights = pop, ad_field = 'ADM1_CODE',
    spatial_obj = NULL, grid = stackers[[z]][[x]])
    tbl$year <- x + 1999
    tbl$mod <- names(stackers)[z]
    return(tbl)
  }))
  return(results)
}))

ggplot() +
  geom_line(data = filter(stack_agg, ADM1_NAME == 'Maharashtra'),
    aes(y = mean, x = year)) +
  facet_wrap(.~mod)
