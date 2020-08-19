library(raster)

grab_ras <- function(indi) {
	setwd("<<<< FILEPATH REDACTED >>>>")
	return(raster(list.files(pattern = '.tif')[2]))
}

indis <- c('w_network', 'w_public_piped', 'w_piped', 'w_imp')
water <- do.call(brick, lapply(indis, grab_ras))
names(water) <- indis
plot(water)

grab_ras <- function(indi) {
	setwd("<<<< FILEPATH REDACTED >>>>")
	file <- list.files(pattern = '.grd')
	print(file)
	return(raster(file, layer = 18))
}
indis <- c('s_network_cr', 's_piped', 's_unimp_cr', 's_imp_cr')
sani <- do.call(brick, lapply(indis, grab_ras))
names(sani) <- indis

spplot(sani)


uncond <- list()
uncond[['s_network']] <- sani[['s_network_cr']] * sani[['s_piped']]
uncond[['s_septic']] <- (1 - sani[['s_network_cr']]) * sani[['s_piped']]
uncond[['s_piped']] <- sani[['s_piped']]
uncond[['s_imp_other']] <- sani[['s_imp_cr']] * (1 - sani[['s_piped']])
uncond[['s_imp_other']] <- sani[['s_imp_cr']] * (1 - sani[['s_piped']])
uncond[['s_imp']] <- sani[['s_imp_cr']] * (1 - sani[['s_piped']]) + (sani[['s_piped']])
uncond[['s_unimp']] <- sani[['s_unimp_cr']] * (1 - uncond[['s_unimp']])
uncond[['s_od']] <- 1 - uncond[['s_imp']] - uncond[['s_unimp']]

spplot(stack(uncond))

setwd("<<<< FILEPATH REDACTED >>>>")
readRDS('subset_shape_dia_south_asia_0')
