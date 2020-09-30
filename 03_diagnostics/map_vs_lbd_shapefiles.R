require(rgdal)
require(ggplot2)

admin_lvl <- c(0:2)

for (admin in admin_lvl){
  map_shp <- readOGR(dsn = '<<< FILEPATH REDACTED >>>', stringsAsFactors = F)
  lbd_shp <- readOGR(dsn = '<<< FILEPATH REDACTED >>>', stringsAsFactors = F)
  pdf('<<< FILEPATH REDACTED >>>', width = 40, height = 20)
  map <- ggplot() + geom_polygon(data = map_shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA, size = 0.1) +
    geom_polygon(data = lbd_shp, aes(x = long, y = lat, group = group), colour = "blue", fill = NA, size = 0.1)
  plot(map)
  dev.off()
}

