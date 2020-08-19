setwd("<<<< FILEPATH REDACTED >>>>")

library(dplyr)
library(raster)

#admin 1 shapefile
ad1 <- shapefile("<<<< FILEPATH REDACTED >>>>")

indis <- c('s_piped', 's_imp_cr', 's_unimp_cr',
  'w_piped', 'w_imp_cr', 'w_unimp_cr')

for (ii in indis) {
	aa <- read.csv(paste0("<<<< FILEPATH REDACTED >>>>"),
	 stringsAsFactors = FALSE)
	aa <- aa[,setdiff(names(aa), names(ad1))]
	d <- SpatialPoints(coords = as.matrix(dplyr::select(aa, longitude, latitude)),
		proj4string = CRS(proj4string(ad1)))
	bb <- sp::over(d, ad1)
	cc <- bind_cols(aa, bb)

	miss_a1 <- filter(cc, is.na(ADM1_CODE))
	nonmiss_a1 <- filter(cc, !is.na(ADM1_CODE))
	for (jj in 1:nrow(miss_a1)) {
		idx <- nonmiss_a1$latitude
		dist <- sqrt((nonmiss_a1$latitude - miss_a1$latitude[jj])^2 + (nonmiss_a1$longitude - miss_a1$longitude[jj])^2)
		latlon <- which(dist == min(dist))
		miss_a1[jj,]$OBJECTID <- nonmiss_a1[latlon,]$OBJECTID
		miss_a1[jj,]$ADM0_NAME <- nonmiss_a1[latlon,]$ADM0_NAME
		miss_a1[jj,]$ADM1_NAME <- nonmiss_a1[latlon,]$ADM1_NAME
		miss_a1[jj,]$ADM1_CODE <- nonmiss_a1[latlon,]$ADM1_CODE
		miss_a1[jj,]$ADM0_CODE <- nonmiss_a1[latlon,]$ADM0_CODE
	}

	cc <- bind_rows(miss_a1, nonmiss_a1)
	print(nrow(filter(cc, is.na(ADM1_CODE))))

	write.csv(cc, file = "<<<< FILEPATH REDACTED >>>>", row.names = FALSE)
}
