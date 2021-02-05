pre_processing <- function(
  rasterPath = "./raster/LC08_L1TP_224079_20201212_20201218_01_T1",
  shapePath = './shp/muni_posadas.shp',
  demPath = "./raster/dem_30m.tif",
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

  }
  
  # Indices espectrales ----
  if (ndvi){
    #For Landsat 8 data, NDVI = (Band 5 – Band 4) / (Band 5 + Band 4)
    ndvi <- ndvi(r_top[['B5_sre']], r_top[['B4_sre']]) 
    names(ndvi) <- 'ndvi'
    r_top <- stack(r_top, ndvi)
  }
  if (evi){
    evi <- overlay(
      r_top[['B5_sre']], r_top[['B4_sre']], 
      r_top[['B2_sre']], 
      fun=function(x,y,z){(2.5*(x-y)/((x+6)*(y-7.5)*(z+1)))})
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
    built <- ndbi(r_top[['B6_sre']], r_top[['B5_sre']])
    names(ndbi) <- 'ndbi'
    r_top <- stack(r_top, ndbi)
  }
  
  #### Corecciones radiometricas TIR
  # Banda TIR de DN a BT
  bt <- br_temp(r$B10_dn, band = "10", conv = TRUE, mult = 0.00033420, add = 0.1, k1 = 774.8853, k2 = 1321.0789)# Tuve que agregar cada uno de los otros parametros
  names(bt) <- 'TIR_bt'
  
  
  # Raster stack de bandas corregidas y TIR 10 en BT
  r_top <- stack(r_top, ndvi, ndwi, built, bt)
  #plotRGB(r_sdos, 8,5,2, stretch = 'lin')
  
  if (!is.null(shapePath)){
    # reading municipality boundary ----
    shape_sf <- read_sf(shapePath)
    
    # croping image ----
    r_sdos_clip <- crop(r_top, shape_sf)
    r_sdos_clip <- mask(r_sdos_clip, shape_sf)
    #plotRGB(r_sdos, 5,4,3, stretch = 'hist')
    writeRaster(r_sdos_clip, paste0(rasterPath, "/", tail(stringr::str_split(rasterPath, '/')[[1]], 1), "_sdos_clip.tif"), overwrite = TRUE)
  }
  
  # saving output ---
  writeRaster(r_sdos, paste0(rasterPath, "/", tail(stringr::str_split(rasterPath, '/')[[1]], 1), "_sdos.tif"), overwrite = TRUE)
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
  mydata <- na.omit(values(r_sdos[[1]]))
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
                  "MacQueen")[1], iter.max = 100)
  
  # Saving statistical results
  save(fit, file = paste0('./outputs/', name , '_KmeansFit.RData'))
  
  # append cluster assignment
  mydata <- data.frame(mydata, fit$cluster) %>% as_tibble()
  
  # creating a raster layer to recieve goup values
  landCover <- r_sdos[[1]]
  # changing pixel values to group values
  landCover[non_na_px] <- mydata$fit.cluster
  
  map <- tm_shape(landCover) +
    tm_raster(palette = 'cat', style = 'cat', legend.show = FALSE) +
    tm_style(style = 'natural')
  tmap_save(map, paste0("./outputs/", name, "_landcover.png"))
  
  writeRaster(landCover, paste0("./raster/", name, "/", name, "_landCover.tif"), overwrite = TRUE)
  
}