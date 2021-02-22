library(sf)
library(ggplot2)
library(tmap)

censo <- read_sf("./shp/Radios Censo 2010/Censo2010Fixed.shp")
censo$area_ha <- st_area(censo)/10000
censo$densidadpobl <- censo$totalpobl/censo$area_ha
plot(st_geometry(censo))
ggplot(censo, aes(x=totalpobl)) + geom_histogram()

pobl <- tm_shape(censo) +
  tm_polygons("totalpobl", style="fisher", title = "Población", border.alpha = 0.1) +
  tm_layout(legend.outside = FALSE, legend.position = c("RIGHT", "BOTTOM")) +
  tm_graticules(lwd = 0) +
  tm_compass(position = c("RIGHT", "TOP"))


densidad <- tm_shape(censo) +
  tm_polygons("densidadpobl", style="fisher", title = "Población/ha", border.alpha = 0.1) +
  tm_layout(legend.outside = FALSE, legend.position = c("RIGHT", "BOTTOM")) +
  tm_graticules(lwd = 0) +
  tm_compass(position = c("RIGHT", "TOP"))

tmap_arrange(pobl, densidad)
sum(censo$totalpobl)
tmap_save(pobl, "./plots/poblacionPosadas.png")
