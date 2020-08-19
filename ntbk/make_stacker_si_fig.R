
stackers <- grab_stacker_raster('w_piped', '2019_05_20_00_00_02',
			'per')



setwd("<<<< FILEPATH REDACTED >>>>")
files <- list.files(pattern = '.grd')
files <- files[grep('_per_0', files)]
inla <- brick(files)

library(ggplot2)
library(rasterVis)
library(viridis)

png("<<<< FILEPATH REDACTED >>>>", units = 'in', 8, 8, res = 300)
print(
gplot(stackers[['gam']][[1]]) +  geom_tile(aes(fill = value)) +
              coord_equal() +
              theme_classic() +
              scale_fill_viridis(begin = 0, end = 1,
                         na.value = 'white',
                         values = c(0, 1), name = 'Access') +
              theme(axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank())
)
dev.off()

png("<<<< FILEPATH REDACTED >>>>", units = 'in', 8, 8, res = 300)
print(

gplot(stackers[['gbm']][[1]]) +  geom_tile(aes(fill = value)) +
              coord_equal() +
              theme_classic() +
              scale_fill_viridis(begin = 0, end = 1,
                         na.value = 'white',
                         values = c(0, 1), name = 'Access') +
              theme(axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank())
)
dev.off()

png("<<<< FILEPATH REDACTED >>>>", units = 'in', 8, 8, res = 300)
print(

gplot(stackers[['lasso']][[1]]) +  geom_tile(aes(fill = value)) +
              coord_equal() +
              theme_classic() +
              scale_fill_viridis(begin = 0, end = 1,
                         na.value = 'white',
                         values = c(0, 1), name = 'Access') +
              theme(axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank())
)
dev.off()

png("<<<< FILEPATH REDACTED >>>>", units = 'in', 8, 8, res = 300)
print(

gplot(inla[[1]]) +  geom_tile(aes(fill = value)) +
              coord_equal() +
              theme_classic() +
              scale_fill_viridis(begin = 0, end = 1,
                         na.value = 'white',
                         values = c(0, 1), name = 'Access') +
              theme(axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                axis.line = element_blank())
)
dev.off()
