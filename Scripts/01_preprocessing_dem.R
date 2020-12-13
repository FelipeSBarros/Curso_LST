library(raster)
library(sf)

# load DEM
dem1 <- raster('./raster/2757-30.img')
dem2 <- raster('./raster/2757-29.img')

# load Posadas
posadas <- read_sf('./shp/muni_posadas.shp')

dem_final <- raster::mosaic(dem1, dem2, fun = mean)

dem_proj <- raster::projectRaster(dem_final, crs = crs(posadas))

plot(dem_proj)
plot(st_geometry(posadas), add = T)

writeRaster(dem_proj, "./raster/dem_30m.tif")
