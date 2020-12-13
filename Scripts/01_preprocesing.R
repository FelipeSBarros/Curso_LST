library(raster)
library(RStoolbox)
library(sf)
library(LSTtools)
library(rgdal)

# loading metadata
m <- list.files(path = "./raster/LC08_L1TP_224079_20200518_20200527_01_T1", pattern ="_MTL.txt$", recursive = TRUE, full.names = TRUE)

m <- readMeta(m)

# loading images ----

# Generar raster stack
r <- stackMeta(m, quantity = "all", category = "image",
               allResolutions = FALSE)

# reading municipality boundary ----
posadas <- read_sf('./shp/muni_posadas.shp')

# croping image ----
r <- crop(r, posadas)
r <- mask(r, posadas)
plotRGB(r, 5,4,3, stretch = 'hist')

# Correcciones ----
h <- estimateHaze(r, hazeBands = 2:6, darkProp = 0.01) # no incluye la banda termica. aca solo calcula la correccion

r_sdos <- radCor(r, metaData = m, method = "sdos", hazeValues = h, bandSet=2:6) # aplica correccion atmosferica; DN a reflectancia
#plotRGB(r_sdos, r=4,g=3,b=2, stretch = "lin")

# Indices espectrales ----
#For Landsat 8 data, NDVI = (Band 5 – Band 4) / (Band 5 + Band 4)
ndvi <- ndvi(r_sdos[['B5_sre']], r_sdos[['B4_sre']]) 
#For Landsat 8 data, NDWI = (Band 5 – Band 6) / (Band 5 + Band 6)
ndwi <- ndvi(r_sdos[['B5_sre']], r_sdos[['B6_sre']]) 
# Calcular built-up
built <- ndbi(r_sdos[['B6_sre']], r_sdos[['B5_sre']])

names(ndvi) <- 'ndvi'
names(ndwi) <- 'ndwi'
names(built) <- 'built'

#### Corecciones radiometricas
# Banda TIR de DN a BT
bt <- br_temp(r$B10_dn, band = "10", conv = TRUE, mult = 0.00033420, add = 0.1, k1 = 774.8853, k2 = 1321.0789)# Tuve que agregar cada uno de los otros parametros
names(bt) <- 'TIR_bt'


# Raster stack de bandas corregidas y TIR 10 en BT
r_sdos <- stack(r_sdos, ndvi, ndwi, built, bt)
#plotRGB(r_sdos, 8,5,2, stretch = 'lin')


# saving output ---
writeRaster(r_sdos, "./raster/LC08_L1TP_224079_20200518_20200527_01_T1_sdos.tif", overwrite = TRUE)





# # Correccion topografica
# dem <- raster("./raster/landsat_xal2/dem.tif")
# summary(dem)
# 
# # Alinear extencion de las imagenes
# dem <- resample(dem, r_clip[[1]], "bilinear")
# 
# r_top <- topCor(r_clip, dem = dem, metaData = m, method = "minnaert") # Corereccion topografica
# plotRGB(r_top,r=4,g=3,b=2,stretch = "lin")



# 
# # Calcular emisividad, obteniendo como output proporcion de la vegetacion
# em <- emissivity(veg, enonveg = 0.95, eveg = 0.99, pveg = TRUE)
# em
# plot(em[[1]], col=colorRampPalette(c("white", "black"))(255))
# plot(em[[2]])
# 
# # Calcular LST
# lst_xal <- landsat_lst(bt, em[[1]], 'L8', conv = FALSE)
# plot(lst_xal, col=brewer.pal(9, 'YlOrRd'))
# plot(ciudad, add=T)
