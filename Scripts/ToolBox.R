pre_processing <- function(
  rasterPath = "./raster/LC08_L1TP_224079_20201212_20201218_01_T1",
  shapePath = './shp/muni_posadas.shp',
  demPath = "./raster/dem_30m.tif",
  ndvi = TRUE,
  evi = TRUE,
  savi = TRUE,
  ndwi = TRUE,
  ndbi = TRUE,
  ...){
  library(raster)
  library(RStoolbox)
  library(sf)
  library(LSTtools)
  library(rgdal)
  
  # loading metadata
  meta <- list.files(path = rasterPath, pattern ="_MTL.txt$", recursive = TRUE, full.names = TRUE)
  
  meta <- readMeta(meta)
  
  # loading images ----
  
  # Generar raster stack
  r <- stackMeta(meta, quantity = "all", category = "image",
                 allResolutions = FALSE)
  
  # Correcciones ----
  h <- estimateHaze(r, hazeBands = 2:6, darkProp = 0.01) # no incluye la banda termica. aca solo calcula la correccion
  
  r_sdos <- radCor(r, metaData = meta, method = "sdos", hazeValues = h, bandSet=2:6) # aplica correccion atmosferica; DN a reflectancia
  #plotRGB(r_sdos, r=4,g=3,b=2, stretch = "lin")
  
  # Correccion topografica 
  # se debe ejecutar 01_preprocessing_dem antes
  if (!is.null(demPath)){
    dem <- raster(demPath)
    #summary(dem)
    
    # Alinear extencion de las imagenes
    dem <- resample(dem, r_sdos[[1]], "bilinear")
    
    r_top <- topCor(r_sdos, dem = dem, metaData = meta, method = "minnaert")

  } else {r_top <- r_sdos}
  
  # Indices espectrales ----
  if (ndvi){
    #For Landsat 8 data, NDVI = (Band 5 – Band 4) / (Band 5 + Band 4)
    ndvi <- ndvi(r_top[['B5_sre']], r_top[['B4_sre']]) 
    names(ndvi) <- 'ndvi'
    r_top <- stack(r_top, ndvi)
  }
  if (evi){
    # In Landsat 8, EVI = 2.5 * ((Band 5 – Band 4) / (Band 5 + 6 * Band 4 – 7.5 * Band 2 + 1)).
    evi <- overlay(
      r_top[['B5_sre']], r_top[['B4_sre']], 
      r_top[['B2_sre']], 
      fun=function(x,y,z){2.5*((x-y)/((x+6*y-7.5*z+1)))})
    names(evi) <- 'evi'
    r_top <- stack(r_top, evi)
  }
  if(savi){
    savi <- overlay(r_top[['B5_sre']], 
                    r_top[['B4_sre']], 
                    fun=function(x,y){((x-y)/(x+y+0.5))*(1.5)})
    names(savi) <- 'savi'
    r_top <- stack(r_top, savi)
  }
  if(ndwi){
    #For Landsat 8 data, NDWI = (Band 5 – Band 6) / (Band 5 + Band 6)
    ndwi <- ndvi(r_top[['B5_sre']], r_top[['B6_sre']]) 
    names(ndwi) <- 'ndwi'
    r_top <- stack(r_top, ndwi)
  }
  if(ndbi){
  # Calcular built-up
    ndbi <- ndbi(r_top[['B6_sre']], r_top[['B5_sre']])
    names(ndbi) <- 'ndbi'
    r_top <- stack(r_top, ndbi)
  }
  
  #### Corecciones radiometricas TIR
  # Banda TIR de DN a BT
  bt <- br_temp(r$B10_dn, band = "10", conv = TRUE, mult = 0.00033420, add = 0.1, k1 = 774.8853, k2 = 1321.0789)# Tuve que agregar cada uno de los otros parametros
  names(bt) <- 'TIR_bt'
  
  
  # Raster stack de bandas corregidas y TIR 10 en BT
  r_top <- stack(r_top, bt)
  #plotRGB(r_sdos, 8,5,2, stretch = 'lin')
  
  if (!is.null(shapePath)){
    # reading municipality boundary ----
    shape_sf <- read_sf(shapePath)
    
    # croping image ----
    r_sdos_clip <- crop(r_top, shape_sf)
    r_sdos_clip <- mask(r_sdos_clip, shape_sf)
    #plotRGB(r_sdos_clip, 5,4,3, stretch = 'hist')
    writeRaster(r_sdos_clip, paste0(rasterPath, "/", tail(stringr::str_split(rasterPath, '/')[[1]], 1), "_sdos_clip.tif"), overwrite = TRUE)
  }
  
  # saving output ---
  writeRaster(r_top, paste0(rasterPath, "/", tail(stringr::str_split(rasterPath, '/')[[1]], 1), "_sdos.tif"), overwrite = TRUE)
}

optimize_groups <- function(initialGroups = 4, 
                            mydata, 
                            seed = 1234, 
                            name){
  library(tibble)
  library(ggplot2)
  
  set.seed(seed)
  wss <- matrix(ncol=3, nrow=2*initialGroups)
  wss[1,1] <- (nrow(mydata) - 1) * sum(apply(mydata,2,var))
  wss[1,2] <- 0
  wss[1,3] <- 0
  
  for (i in 2:nrow(wss)){
    set.seed(seed)
    wss[i, 1] <- sum(kmeans(mydata, centers = i)$withinss)
    wss[i, 2] <- abs(wss[i-1, 1] - wss[i, 1])
    wss[i, 3] <- wss[i, 2] / wss[2, 2]
    #wss[i, 4] <- wss[i, 1] / wss[1, 1]
  }
  
  wss <- as_tibble(wss)
  names(wss) <- c('wss', 'DeltaChange', 'RatioChange')
  wss$group <- as.integer(rownames(wss))
  
  plot <- ggplot(wss, aes(x = as.factor(group), y = wss, group = 1)) + 
    geom_line() + 
    labs(x = "Number of Clusters", 
         y = "Within groups sum of squares") + 
    theme(text = element_text(size = 17)) + 
    geom_point() + 
    annotate("rect",
             xmin = which(wss$RatioChange==min(wss$RatioChange[-1]))-1,
             xmax = which(wss$RatioChange==min(wss$RatioChange[-1])),
             ymin = min(wss$wss), 
             ymax = wss[which(wss$RatioChange==min(wss$RatioChange[-1])),'wss'][[1]], 
             alpha=.3, 
             fill = "red") +
    annotate("text", 
             x = which(wss$RatioChange==min(wss$RatioChange[-1]))-1, 
             y = wss[which(wss$RatioChange==min(wss$RatioChange[-1])),'wss'][[1]], 
             label = paste(which(wss$RatioChange==min(wss$RatioChange[-1]))-1), 
             vjust=-1.5, 
             size=5, 
             colour='red')
  
  
  if (!file.exists("./plots/")) dir.create("./plots/")
  ggsave(paste0("./plots/", name, "_Kmeans_clusterAnalysis.png"), dpi = 300)
  # dev.off()
  }

landcover <- function(
  rasterPath = "./raster/LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip.tif",
  shapePath = './shp/muni_posadas.shp',
  initialGroups = 4
  ){
  library(sf)
  library(raster)
  library(tibble)
  library(tmap)
  library(dplyr)
  library(tidyr)
  
  # Getting project name
  name <- sub('.tif','',
              tail(stringr::str_split(
                rasterPath, '/')[[1]], 1))
  
  # creating output folder
  if (!file.exists("./outputs")){
    dir.create("./outputs")}
  
  # loading data ----
  r_sdos <- stack(rasterPath)
  
  shape_sf <- read_sf(shapePath)
  
  # Processing data ----
  # Identifying pxls w/ NA
  non_na_px <- which(!is.na(getValues(r_sdos[[1]])))
  
  # Creating data to Kmeans analysis
  mydata <- na.omit(values(r_sdos))
  mydata <- scale(mydata) %>% as_tibble() # standardize variables

  # optimizing amoutn of groups ----
  optimize_groups(mydata = mydata, name = name, initialGroups)
  nClasses <- as.integer(readline("Define amount of classes t be created:"))
  # training algorithm for unsupervised ----
  # K-Means Cluster Analysis
  fit <- kmeans(
    mydata, nClasses, 
    algorithm = c("Hartigan-Wong",
                  "Lloyd",
                  "Forgy",
                  "MacQueen")[1], iter.max = 300)
  
  # Saving statistical results
  save(fit, file = paste0('./outputs/', name , '_KmeansFit.RData'))
  
  # append cluster assignment
  mydata <- data.frame(mydata, fit$cluster) %>% as_tibble()
  
  # Creating spectral sign plot
  head(mydata)
  names(mydata) <- c(paste0("B", 1:nlayers(r_sdos)), "Class")
  t <- tidyr::pivot_longer(mydata, !Class, names_to = "bands", values_to = "RT") %>% group_by(bands, Class) %>% 
    summarise(
      mean = mean(RT),
      sd = sd(RT),
      # n = sum(RT),
      # se = sd/sqrt(n)
    ) #%>% filter(bands != "B10")
  
  spectralSign <- ggplot(t, aes(x=as.factor(bands),
                                y=mean, 
                                group = as.factor(Class))) + 
    geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), alpha=0.2) + 
    geom_line(
      aes(colour = as.factor(Class))) +
    labs(x = "Satelitte image bands and spectral indices", 
         y = "Reflectance") + 
    theme(text = element_text(size = 17)) + guides(colour=guide_legend(title=NULL))
  if (!file.exists("./plots/")) dir.create("./plots/")
  ggsave(paste0("./plots/", name, "_SpectralSign.png"), dpi = 300)
  
  
  density <- mydata %>% tidyr::pivot_longer(!Class, names_to = "bands", values_to = "RT") %>% 
    ggplot(., aes(x=RT, group=as.factor(Class), fill=as.factor(Class))) +
    geom_density(alpha = 0.75) + 
    xlim(-2.75, 2.75) +
    geom_vline(
      data = . %>% group_by(Class) %>% summarise(grp.mean = mean(RT)),
      aes(xintercept=grp.mean, color = as.factor(Class)), linetype="dashed", size=1) + theme(panel.background = element_blank(),
          panel.grid.major = element_line(color = "gray", size = 0.5),
          panel.grid.minor = element_line(color = "gray", size = 0.5),
          axis.ticks = element_blank()) +
    labs(x = "Reflectance Value",
         y = "Density",
         title = "Density histograms of spectral profiles",
         subtitle = "Vertical lines represent mean group reflectance values")
  ggsave(paste0("./plots/", name, "_SpectralDensity.png"), dpi = 300)
  
  
  # creating a raster layer to recieve goup values
  landCover <- r_sdos[[1]]
  # changing pixel values to group values
  landCover[non_na_px] <- mydata$Class
  
  map <- tm_shape(landCover) +
    tm_raster(palette = 'cat', style = 'cat', legend.show = TRUE, title = "Landcover class",
              labels = c("Foresta", 
                         "Suelo Expuesto",
                         "Pastizales",
                         "Agua",
                         "Área Urbana")) +
    tm_layout(legend.outside = FALSE, 
              legend.position = c("RIGHT", "BOTTOM")) +
    tm_graticules(lwd = 0) #+
    tm_compass(position = c("RIGHT", "TOP"))
  
  tmap_save(map, paste0("./plots/", name, "_landcover.png"))
  
  rgb <- r_sdos[[1:8]] %>% stretch()
  rgb <- tm_shape(rgb) +
    tm_rgb(r = 4, g=3, b=1) +
    tm_graticules(lwd = 0) +
    tm_compass()
  
  composition <- tmap_arrange(rgb, map)
  tmap_save(composition, paste0("./plots/", name, "_composition.png"))
  
  strsplited <- stringr::str_split(name, "_")
  namePath <- paste(strsplited[[1]][1:(length(strsplited[[1]])-2)], collapse = '_')
  
  writeRaster(landCover, paste0("./raster/", namePath, "/", name, "_landCover.tif"), overwrite = TRUE)
  # landCover <- raster('./raster/LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip_landCover.tif')
}

automate_suhi <- function(
  rasterPath = "./raster/LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip.tif",
  ndviLayer = 6,
  ndbiLayer = 10, 
  btLayer = 11, # band temperature layer
  urbanID = 5,
  shapePath = './shp/muni_posadas.shp',
  landcoverPath = './raster/LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip_landCover.tif'){
  library(sf)
  library(raster)
  library(LSTtools)
  library(tmap)
  library(ggplot2)
  library(usdm)
  library(tibble)
  library(readr)
  
  # defining project name
  name <- sub('.tif','',
              tail(stringr::str_split(
                rasterPath, '/')[[1]], 1))
  strsplited <- stringr::str_split(name, "_")
  namePath <- paste(strsplited[[1]][1:(length(strsplited[[1]])-2)], collapse = '_')
  
  # loading data ----
  r_sdos <- stack(rasterPath)
  posadas <- read_sf(shapePath)
  st <- raster::stretch(r_sdos)
  rgb <- r_sdos[[1:8]] %>% stretch()

  # leer land cover
  cover_r <- raster(landcoverPath)
  
  # emissivity calculation based on NDVI----
  em <- emissivity(r_sdos[[ndviLayer]], pveg = TRUE)
  
  # maps
  emissivity <- tm_shape(em[[1]]) +
    tm_raster(style='fisher') +
    tm_layout(legend.outside = FALSE, legend.position = c("RIGHT", "BOTTOM")) +
    tm_graticules(lwd = 0)
  #tmap_save(emissivity, paste0("./plots/", name, "_emissivity.png"))
  
  pveg <- tm_shape(em[[2]]) +
    tm_raster(style='fisher', title = "Proportion of vegetation") +
    tm_layout(legend.outside = FALSE, legend.position = c("RIGHT", "BOTTOM")) +
    tm_graticules(lwd = 0)
  #tmap_save(pveg, paste0("./plots/", name, "_pveg.png")
  
  both <- tmap_arrange(emissivity, pveg)
  tmap_save(both, paste0("./plots/", name, "_em_pveg.png"))
  
  # Calcular <- ----
  lst <- landsat_lst(r_sdos[[btLayer]], em[[1]], 'L8', conv = FALSE)
  
  # saving raster result
  writeRaster(lst, paste0("./raster/", namePath, "/", name, "_lst.tif"), overwrite = TRUE)
  # lst <- raster("./raster/LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip_lst.tif")
  
  # saving map
  lstMap <- tm_shape(lst) + 
    tm_raster(style = 'fisher', title = 'Temperatura °C', palette = 'YlOrRd') +
    tm_layout(legend.outside = FALSE, 
              legend.position = c("RIGHT", "BOTTOM"),
              main.title = "Land Surface Temperature") +
    tm_graticules(lwd = 0)
  tmap_save(lstMap, paste0("./plots/", name, "_LST.png"))
  
  #Calcular estadisticas zonales ----
  
  # Urban Heat Island UHI ----
  # Calcular estadisticas basicas
  st <- uhi_stats(lst, cover_r, id=urbanID) #calcula las diferencias con el id dos, area urbana
  st$clase <- c('Floresta', 'Suelo expuesto', 'Pastizales', 'Agua', 'Area Urbana')
  write.csv(st, paste0("./outputs/", name, "_uhi_stats.csv"), row.names = FALSE)
  
  # Calcular hot island area HIA ----
  if (!is.null(urbanID)){
    print("HIA for Urban ID")
    lst_city <- lst * (cover_r == urbanID) # duvida: HIa so para ciudad
    lst_city[lst_city==0] <- NA
  } else {
    print("HIA for the entire layer")
    lst_city <- lst
  }
  
  #plot(lst_city)
  # HIA
  #hia_cal <- raster("./raster/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip_hia.tif")
  hia_cal <- hia(lst_city)
  write_rds(hia_cal, paste0("./outputs/", name, "_hia"))
  
  # saving map
  hia <- tm_shape(rgb) +
    tm_rgb(4,3, 2, alpha = 0.5) +
    tm_shape(hia_cal) + 
    tm_raster(style = 'fisher', title = 'Temperatura °C', palette = 'YlOrRd') +
    tm_layout(legend.outside = FALSE, 
              legend.position = c("RIGHT", "BOTTOM"),
              main.title = "Islas de calor urbano - Posadas, Misiones",
              main.title.size = .9) +
    tm_graticules(lwd = 0)
  tmap_save(hia, paste0('./plots/', name, '_hiaMap.png'))
  writeRaster(hia_cal[[1]], paste0('./raster/', name, '_hia.tif'), overwrite = TRUE)
  # hia_cal <- raster("./raster/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip_hia.tif")
  
  # Correlacion ----

  # spots <- getis(lst, dist = 70, p=0.05)
  print("running LISA")
  g <- lisa(lst, d1 = 0, d2 = 70, statistic = "G*")
  print("LISA done")
  Sys.sleep(10)
  # tm_shape(g) +
  #   tm_raster(style = "fisher")
  
  print("Converting to polygon")
  pv <- round(2 * pnorm(-abs(getValues(g))), 7)
  pv <- as.data.frame(pv)
  FDR <- round(p.adjust(pv[, 1], "fdr"), 7)
  g <- round(g, 2)
  Sys.sleep(10)
  writeRaster(g, paste0("./raster/", namePath, "/", name, "_lisa.tif"), overwrite = TRUE)
  # g <- raster("./raster/LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip_lisa.tif")
  Sys.sleep(10)
  # gdal_polygonize("./raster/LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip_lisa.tif", "./shp//LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip_lisa.shp")
  #v.g <- clamp(g)
  v.g <- rasterToPolygons(g, na.rm = TRUE)
  v.g <- st_as_sf(v.g)
  names(v.g)[1] <- "Z.scores"
  v.g$FDR <- na.omit(FDR)
  v.g$cluster <- ifelse(v.g$Z.scores > 0 & v.g$FDR < 0.05, "Hot spot", ifelse(v.g$Z.scores < 0 & v.g$FDR < 0.05, "Cold spot", "No sig"))
  
  write_sf(st_as_sf(v.g), paste0("./shp/", name, "_getis_analysis.shp"), overwrite = TRUE)
  # v.g <- read_sf("./shp/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip_getis_analysis.shp")
  Sys.sleep(10)
  
  print("Mapa Area correlação")
  AreasCorrelacion <- tm_shape(v.g) +
    tm_fill("cluster") +
    tm_layout(legend.outside = FALSE, 
            legend.position = c("RIGHT", "BOTTOM"),
            main.title = "D",
            main.title.size = .9) +
    tm_graticules(lwd = 0)
  tmap_save(AreasCorrelacion, paste0("./plots/", name, "_AreasCorrelacionG.png"))
  
  indiceG <- tm_shape(g) +
    tm_raster(style = 'fisher') +
    tm_layout(legend.outside = FALSE,
              legend.position = c("RIGHT", "BOTTOM"),
              main.title = "Indice G",
              main.title.size = .9) +
    tm_graticules(lwd = 0)
  tmap_save(indiceG, paste0("./plots/", name, "_indiceG.png"))
  
  
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
    st_indices <- stack(lst, r_sdos[[ndviLayer]], r_sdos[[ndbiLayer]])
  names(st_indices) <- c("LST", "NDVI", "NDBI")
  # st_indices <- mask(st_indices, ciudad)
  
  
  table <- as.data.frame(na.omit(values(st_indices)))
  table <- as_tibble(table)
  
  
  # ggplot(table, aes(x = NDVI, y = LST)) + 
  #   geom_point(color='grey') + 
  #   geom_smooth(method = "lm", color = 'red') + theme_minimal() +
  #   labs(title = "Grafico dispersión NDVI~LST")
  # ggsave(paste0('./plots/', name, 'NDVI_LST.png')
  
  # ggplot(table, aes(x = NDVI, y = LST)) + 
  #   geom_point(color='grey') + 
  #   geom_smooth(color = 'red') + theme_minimal() +
  #   labs(title = "Grafico dispersión NDVI~LST")
  #ggsave(paste0('./plots/', name, '_NDVI_LST_gam.png'))
  
  
  # ggplot(table, aes(x = NDBI, y = LST)) + 
  #   geom_point(color='lightgrey', alpha = .5) + 
  #   geom_smooth(color = 'red') + theme_minimal() +
  #   labs(title = "Grafico dispersión NDBI~LST")
  # ggsave(paste0('./plots/', name, '_NDBI_LST.png'))
  
  
  # Estructura de los valores
  # ggplot(table, aes(NDVI)) + geom_histogram() + theme_minimal() +
  #   labs(title = "Histograma de NDVI")
  # ggplot(table, aes(NDBI)) + geom_histogram() + theme_minimal() +
  #   labs(title = "Histograma de NDBI")

  
  # Correlacion
  correlacion <- layerStats(st_indices, 'pearson', na.rm=T)
  cor <- as_tibble(correlacion$'pearson correlation coefficient')
  write_csv(cor, paste0('./outputs/', name, '_correlacion.csv'))
  
  cormean <- as.data.frame(correlacion$mean)
  cormean$index <- rownames(cormean)
  names(cormean)[1] <- 'mean'
  write_csv(cormean, paste0('./outputs/', name, '_correlacion_mean.csv'))
}

landcover2vec <- function(landcoverPath = "./raster/LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip_landCover.tif"){
  library(raster)
  library(sf)
  
  # Getting project name
  name <- sub('.tif','',
              tail(stringr::str_split(
                landcoverPath, '/')[[1]], 1))
  landcover.r <- raster(landcoverPath)
  names(landcover.r) <- 'classes'
  landcover.v <- rasterToPolygons(landcover.r)
  landcover.v <- st_as_sf(landcover.v)
  landcover.v <- landcover.v %>% 
    dplyr::group_by(classes) %>%
    dplyr::summarize()
  landcover.v <- st_cast(landcover.v,"POLYGON")
  
  write_sf(landcover.v, paste0("./shp/", name, ".shp"))
}

lst_explore <- function(lstPath = "./raster/LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip_lst.tif", rasterPath = "./raster/LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip.tif",
               landcoverPath = './raster/LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip_landCover.tif',
               specIndexes = c(6, 7, 8, 9, 10)){
  # HIA
  r <- raster(lstPath)
  names(r) <- "LST"
  # landcover
  landcover <- raster(landcoverPath)
  #landcover <- landcover * (landcover == 5)
  #landcover[landcover==0] <- NA
  names(landcover) <- "landcover"
  
  # spectral indexes
  stackRaster <- stack(rasterPath)[[specIndexes]]
  r <- projectRaster(r, stackRaster[[1]])
  #stackRaster <- stackRaster * (landcover == 5)
  names(stackRaster) <- c("NDVI", "EVI","SAVI", "NDWI", "NDBI")
  # compareRaster(r, landcover, extent = T, crs = T, res = T, orig = T)
  
  r <- addLayer(r, stackRaster, landcover)
  r.df <- na.omit(as.data.frame(r))
  #head(r.df)
  write.csv(r.df, "./outputs/Hia_specIndexes.csv", row.names = FALSE)
}
