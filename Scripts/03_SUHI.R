library(sf)
library(raster)
library(LSTtools)
library(tmap)
library(ggplot2)
library(usdm)
library(tibble)

# loading data ----
r_sdos <- stack("./raster/LC08_L1TP_224079_20200518_20200527_01_T1_sdos.tif")
posadas <- read_sf('./shp/muni_posadas.shp')
# leer land cover
cover_r <- raster('./raster/landCover.tif')

# emissivity calculation ----
em <- emissivity(r_sdos[[6]], pveg = TRUE)

# maps
emissivity <- tm_shape(em[[1]]) +
  tm_raster(style='fisher') +
  tm_layout(main.title = "Emissivity", legend.outside = TRUE) +
  tm_graticules(lwd = 0)
#tmap_save(emissivity, "./outputs/emissivity.png")

pveg <- tm_shape(em[[2]]) +
  tm_raster(style='fisher') +
  tm_layout(main.title = "Proportion of vegetation", legend.outside = TRUE) +
  tm_graticules(lwd = 0)
#tmap_save(pveg, "./outputs/pveg.png")

both <- tmap_arrange(emissivity, pveg)
tmap_save(both, "./outputs/em_pveg.png")

# Calcular LST ----
lst_posadas <- landsat_lst(r_sdos[[9]], em[[1]], 'L8', conv = FALSE)

# saving raster result
#writeRaster(lst_posadas, "./raster/lst_posadas.tif", overwrite = TRUE)

# saving map
lst <- tm_shape(lst_posadas) + 
  tm_raster(style = 'fisher', title = 'Temperatura *C', palette = 'YlOrRd') +
  tm_layout(legend.outside=T, main.title = "Land Surface Temperature - Posadas, Misiones") +
  tm_graticules(lwd = 0)
lst
#tmap_save(lst, './outputs/LST_posadas.png')

#Calcular estadisticas zonales ----

# Urban Heat Island UHI ----
# Calcular estadisticas basicas
st <- uhi_stats(lst_posadas, cover_r, id=5) #calcula las diferencias con el id dos, area urbana
st$clase <- c('Herbacea/Pasto', 'Agua', 'Foresta', 'Suelo exposto', 'Area urbana')
#write.csv(st, "./outputs/uhi_stats.csv", row.names = FALSE)

# Calcular hot island area HIA ----
plot(cover_r == 5)
lst_city <- lst_posadas# * (cover_r == 5) # duvida: HIa so para ciudad
plot(lst_city)
# HIA
hia_cal <- hia(lst_city)
#write_csv(as_tibble(hia_cal[2:length(hia_cal)]), "./outputs/hia.csv")

# saving map
hia <- tm_shape(hia_cal[[1]]) + 
  tm_raster(style = 'fisher', title = 'Temperatura *C', palette = 'YlOrRd') +
  tm_layout(legend.outside=T, main.title = "Areas de Islas de calor urbana") +
  tm_graticules(lwd = 0)
hia
# tmap_save(hia, './outputs/HIA_posadas_AU.png')

st <- raster::stretch(r_sdos)
hiaMap <- tm_shape(st) + 
  tm_rgb(5,4,3) +
  tm_shape(hia_cal[[1]]) + 
  tm_raster(style = 'fisher', title = 'Temperatura *C', palette = 'YlOrRd', colorNA = NULL) +
  tm_layout(legend.outside=T, main.title = "Areas de Islas de calor urbana") +
  tm_graticules(lwd = 0)

#tmap_save(hiaMap, './outputs/hiaMap.png')

# Correlacion ----

#spots <- getis(lst_xal, dist = 70, p=0.05) 
#g <- lisa(lst_xal, d1 = 0, d2 = 70, statistic = "G*")
#pv <- round(2 * pnorm(-abs(getValues(g))), 7)
#pv <- as.data.frame(pv)
#FDR <- round(p.adjust(pv[, 1], "fdr"), 7)
#g <- round(g, 2)
#writeRaster(g, "./raster/g.tif")
# v.g <- clamp(g)
# v.g <- rasterToPolygons(v.g, na.rm = TRUE)

#names(v.g) <- "Z.scores"
#v.g$FDR <- na.omit(FDR)
#v.g$cluster <- ifelse(v.g$Z.scores > 0 & v.g$FDR < 0.05, "Hot spot", ifelse(v.g$Z.scores < 0 & v.g$FDR < 0.05, "Cold spot", "No sig"))

#write_sf(st_as_sf(v.g), "./shp/getis_analysis.shp", overwrite = TRUE)

g <- raster("./raster/g.tif")
v.g <- read_sf("./shp/getis_analysis.shp")

tm_shape(v.g) +
  tm_fill("cluster")

indiceG <- tm_shape(g)+
  tm_raster(style = 'fisher') +
  tm_layout(legend.outside=T, main.title = "Indice G de correlación espacial de temperatura") 
#tmap_save(indiceG, "./outputs/indiceG.png")

# Aggregate by cluster type and plot with proper colors
# spots_a <- sf::aggregate(v.g, by = "cluster")
# spots_a$color <- ifelse(spots_a$cluster == "Hot spot", "red", 
#                         ifelse(spots_a$cluster == "Cold spot", "blue", "grey"))
# plot(spots_a, col=spots_a$color, border=NA, axes=F, main = "Hot-cold spots")  
#   legend("bottom", ncol = 3, fill=spots_a$color, legend = c("Cold", "Hot", "No sig."))
# png("./outputs/result.png")
# dev.off()

# Calcular SUHI ----
# Extraer solo buffer
# plot(ciudad)
# plot(b, col = "black")
# buf <- b - ciudad
# plot(buf, col = "black")

# Nivel ciudad, todas las clases
# xal_mean <- cellStats(mask(lst_xal, ciudad), 'mean')
# buf_mean <- cellStats(mask(lst_xal, buf), 'mean')
# suhi <- # Formula SUHI
#   suhi

# Solo pixeles impermeables
# imp_mean <- cellStats(mask(lst_posadas, cover[cover$ID == 5,]), 'mean')
# suhi <- # Formula SUHI
#   suhi

# Cada pixel de la ciudad vs referencia
# suhi_pix <- mask(lst_xal, ciudad) - buf_mean
# plot(suhi_pix, col=brewer.pal(9, 'YlOrRd'))

# Controlando elevacion
# Identificar rango de elevacion en la ciudad
# dem <- raster("./raster/dem_30m.tif")
# dem_ciudad <- trim(mask(dem, ciudad))
# plot(dem_ciudad)
# elev_ciudad <- range(dem_ciudad[], na.rm=TRUE)
# elev_ciudad

# Rango de valores del DEM en el buffer
# range(mask(dem, buf)[], na.rm=TRUE)

# Crear capa de DEM del buffer
# dem_buf <- mask(dem, buf)
# plot(dem_buf)

# eliminar pixeles de acuerdo a criterio de exclusion
# dem_buf[dem_buf < min(elev_ciudad) - 50 | dem_buf > max(elev_ciudad) + 50] <- NA
# plot(dem_buf)

# Calcular SUHI controlando elevacion
# buf_mean <- cellStats(mask(lst_xal, dem_buf), 'mean')
# suhi_elev <- # Formula SUHI
#   suhi_elev

# Plot NDVI vs LST
st_indices <- stack(lst_posadas, r_sdos[[6]], r_sdos[[7]], r_sdos[[8]])
names(st_indices) <- c("LST", "NDVI", "NDWI", "NDBI")
# st_indices <- mask(st_indices, ciudad)


table <- as.data.frame(na.omit(values(st_indices)))
table <- as_tibble(table)


ggplot(table, aes(x = NDVI, y = LST)) + 
  geom_point(color='grey') + 
  geom_smooth(method = "lm", color = 'red') + theme_minimal() +
  labs(title = "Grafico dispersión NDVI~LST")
ggsave('./outputs/NDVI_LST.png')

ggplot(table, aes(x = NDVI, y = LST)) + 
  geom_point(color='grey') + 
  geom_smooth(color = 'red') + theme_minimal() +
  labs(title = "Grafico dispersión NDVI~LST")
ggsave('./outputs/NDVI_LST_gam.png')


ggplot(table, aes(x = NDBI, y = LST)) + 
  geom_point(color='lightgrey', alpha = .5) + 
  geom_smooth(method = "lm", color = 'red') + theme_minimal() +
  labs(title = "Grafico dispersión NDBI~LST")
ggsave('./outputs/NDBI_LST.png')

ggplot(table, aes(x = NDWI, y = LST)) + 
  geom_point(color='lightgrey', alpha = .5) + 
  geom_smooth(method = "lm", color = 'red') + theme_minimal() +
  labs(title = "Grafico dispersión NDWI~LST")
ggsave('./outputs/NDWI_LST.png')

# Estructura de los valores
ggplot(table, aes(NDVI)) + geom_histogram() + theme_minimal() +
  labs(title = "Histograma de NDVI")
ggplot(table, aes(NDBI)) + geom_histogram() + theme_minimal() +
  labs(title = "Histograma de NDBI")
ggplot(table, aes(NDWI)) + geom_histogram() + theme_minimal() +
  labs(title = "Histograma de NDWI")

# Correlacion
correlacion <- layerStats(st_indices, 'pearson', na.rm=T)
cor <- as_tibble(correlacion$'pearson correlation coefficient')
write_csv(cor, "./outputs/correlacion.csv")

cormean <- as.data.frame(correlacion$mean)
cormean$index <- rownames(cormean)
names(cormean)[1] <- 'mean'
write_csv(cormean, "./outputs/correlacion_mean.csv")
