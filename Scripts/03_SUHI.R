library(sf)
library(raster)
library(LSTtools)
library(tmap)
# library(RStoolbox)
# library(RColorBrewer)
# library(rgdal)
library(ggplot2)
library(usdm)

# loading data ----
r_sdos <- stack("./raster/LC08_L1TP_224079_20200518_20200527_01_T1_sdos.tif")
posadas <- read_sf('./shp/muni_posadas.shp')
# leer land cover
cover_r <- raster('./raster/landCover.tif')

# emissivity calculation ----
em <- emissivity(r_sdos[[6]], veg = 0.8, nonveg = 0.7, enonveg = 0.7, eveg = 0.8, pveg = TRUE)
plot(em[[1]], col=colorRampPalette(c("white", "black"))(255), main = 'Emissivity')
plot(em[[2]], main = 'Proportion of vegetation')

# Calcular LST ----
lst_posadas <- landsat_lst(r_sdos[[9]], em[[1]], 'L8', conv = FALSE)
# saving raster result
writeRaster(lst_posadas, "./raster/lst_posadas.tif")
# saving map
lst <- tm_shape(lst_posadas) + 
  tm_raster(style = 'fisher', title = 'Temperatura *C', palette = 'YlOrRd') +
  tm_layout(legend.outside=T, main.title = "Land Surface Temperature - Posadas, Misiones") +
  tm_graticules(lwd = 0)
lst
tmap_save(lst, './outputs/LST_posadas.png')

#Calcular estadisticas zonales ----

# Urban Heat Island UHI ----
# Calcular estadisticas basicas
st <- uhi_stats(lst_posadas, cover_r, id=5) #calcula las diferencias con el id dos, area urbana
st$clase <- c('Herbacea/Pasto', 'Agua', 'Foresta', 'Suelo exposto', 'Area urbana')
write.csv(st, "./outputs/uhi_stats.csv", row.names = FALSE)

# Calcular hot island area HIA ----
plot(cover_r == 5)
lst_city <- lst_posadas * (cover_r == 5) # duvida: HIa so para ciudad
plot(lst_city)
hia_cal <- hia(lst_city)
hia_cal
# saving map
hia <- tm_shape(hia_cal[[1]]) + 
  tm_raster(style = 'fisher', title = 'Temperatura *C', palette = 'YlOrRd') +
  tm_layout(legend.outside=T, main.title = "Areas de Islas de calor urbana") +
  tm_graticules(lwd = 0)
hia
tmap_save(hia, './outputs/HIA_posadas_AU.png')

# Correlacion ----
plot(lst_posadas)
#spots <- getis(lst_xal, dist = 70, p=0.05) 
#g <- lisa(lst_xal, d1 = 0, d2 = 70, statistic = "G*")
#pv <- round(2 * pnorm(-abs(getValues(g))), 7)
#pv <- as.data.frame(pv)
#FDR <- round(p.adjust(pv[, 1], "fdr"), 7)
#g <- round(g, 2)
#writeRaster(g, "./raster/g.tif")

#v.g <- rasterToPolygons(g, na.rm = TRUE)

#names(v.g) <- "Z.scores"
#v.g$FDR <- na.omit(FDR)
#v.g$cluster <- ifelse(v.g$Z.scores > 0 & v.g$FDR < 0.05, "Hot spot", ifelse(v.g$Z.scores < 0 & v.g$FDR < 0.05, "Cold spot", "No sig"))

#write_sf(st_as_sf(v.g), "./shp/getis_analysis.shp")

g <- raster("./raster/g.tif")
v.g <- read_sf("./shp/getis_analysis.shp")

# Aggregate by cluster type and plot with proper colors
spots_a <- sf::aggregate(v.g, by = "cluster")
spots_a$color <- ifelse(spots_a$cluster == "Hot spot", "red", 
                        ifelse(spots_a$cluster == "Cold spot", "blue", "grey"))
plot(spots_a, col=spots_a$color, border=NA, axes=F, main = "Hot-cold spots")  
  legend("bottom", ncol = 3, fill=spots_a$color, legend = c("Cold", "Hot", "No sig."))
png("./outputs/result.png")
dev.off()

# Calcular SUHI ----
# Extraer solo buffer
plot(ciudad)
plot(b, col = "black")
buf <- b - ciudad
plot(buf, col = "black")

# Nivel ciudad, todas las clases
xal_mean <- cellStats(mask(lst_xal, ciudad), 'mean')
buf_mean <- cellStats(mask(lst_xal, buf), 'mean')
suhi <- # Formula SUHI
  suhi

# Solo pixeles impermeables
imp_mean <- cellStats(mask(lst_xal, cover[cover$ID == 2,]), 'mean')
suhi <- # Formula SUHI
  suhi

# Cada pixel de la ciudad vs referencia
suhi_pix <- mask(lst_xal, ciudad) - buf_mean
plot(suhi_pix, col=brewer.pal(9, 'YlOrRd'))

# Controlando elevacion
# Identificar rango de elevacion en la ciudad
dem_ciudad <- trim(mask(dem, ciudad))
plot(dem_ciudad)
elev_ciudad <- range(dem_ciudad[], na.rm=TRUE)
elev_ciudad

# Rango de valores del DEM en el buffer
range(mask(dem, buf)[], na.rm=TRUE)

# Crear capa de DEM del buffer
dem_buf <- mask(dem, buf)
plot(dem_buf)

# eliminar pixeles de acuerdo a criterio de exclusion
dem_buf[dem_buf < min(elev_ciudad) - 50 | dem_buf > max(elev_ciudad) + 50] <- NA
plot(dem_buf)

# Calcular SUHI controlando elevacion
buf_mean <- cellStats(mask(lst_xal, dem_buf), 'mean')
suhi_elev <- # Formula SUHI
  suhi_elev

# Plot NDVI vs LST
st_indices <- stack(lst_xal, veg, built)
names(st_indices) <- c("LST", "NDVI", "NDBI")
st_indices <- mask(st_indices, ciudad)
scatter.smooth(x=values(st_indices[[2]]), y=values(st_indices[[1]]), ylab = "LST  (?C)",
               main="LST ~ NDVI", col = alpha("grey23", 0.15), xlab = "NDVI")
scatter.smooth(x=values(st_indices[[3]]), y=values(st_indices[[1]]), ylab = "LST  (?C)",
               main="LST ~ NDBI", col = alpha("grey23", 0.15), xlab = "NDBI")

# Estructura de los valores
hist(st_indices[[1]])
hist(st_indices[[2]])

# Correlacion
correlacion <- layerStats(st_indices, 'pearson', na.rm=T)
corr_matrix <- correlacion$'pearson correlation coefficient'
corr_matrix
