rm(list = ls())
library(raster)

indicator <- 'w_imp_cr'
reg <- 'dia_nessa'
rd <- '2020_01_01_00_00_82'

setwd("<<<< FILEPATH REDACTED >>>>")
world_map <- brick(paste0(indicator, '_prediction_eb_bin0_', reg, '_0.grd'))

setwd("<<<< FILEPATH REDACTED >>>>")
load(paste0(rd, '_bin0_', reg, '_0.RData'))

spplot(stack(world_map[[18]], cov_list$gam[[18]],
  cov_list$gbm[[18]], cov_list$lasso[[18]]))

library(data.table)
library(dplyr)
library(ggplot2)
library(rasterVis)
library(RColorBrewer)
library(ggpubr)
mydat <- fread("<<<< FILEPATH REDACTED >>>>")

tza <- filter(mydat, country %in% c('KEN', 'UGA', 'TZA'))

gglist <- list()

gglist[[1]] <- gplot(world_map[[18]]) +
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = 'Greys') +
  scale_color_distiller(palette = 'RdYlBu') +
  ggtitle('INLA')

gglist[[2]] <- gplot(cov_list$gam[[18]]) +
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = 'Greys') +
  scale_color_distiller(palette = 'RdYlBu') +
  ggtitle('GAM')

gglist[[3]] <- gplot(cov_list$gbm[[18]]) +
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = 'Greys') +
  scale_color_distiller(palette = 'RdYlBu') +
  ggtitle('GBM')


gglist[[4]] <- gplot(cov_list$lasso[[18]]) +
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = 'Greys') +
  scale_color_distiller(palette = 'RdYlBu') +
  ggtitle('LASSO')

ggarrange(plotlist = gglist)

gplot(world_map[[18]]) +
  geom_tile(aes(fill = value)) +
  geom_point(data = tza, aes(x = longitude, y = latitude,
    col = w_piped/N),
    alpha = 0.5) +
  scale_fill_distiller(palette = 'Greys') +
  scale_color_distiller(palette = 'RdYlBu') +
  ggtitle('INLA')


gglist[[1]]



setwd("<<<< FILEPATH REDACTED >>>>")
list.files(pattern = 'inla')
r <- readRDS("<<<< FILEPATH REDACTED >>>>")
pgrid0 <- inla.mesh.projector(mesh_s, dims=c(2, 1800))
prd0 <- inla.mesh.project(pgrid0, r$summary.ran$sp.no.t$mean)

prd_test <- prd0[which(apply(prd0, 1, sum, na.rm = TRUE) != 0),which(apply(prd0, 2, sum, na.rm = TRUE) != 0)]

levelplot(prd_test, col.regions=topo.colors(99),
  main='latent field mean: binomial', xlab='', ylab='', scales=list(draw=FALSE))

levelplot(prd0)

