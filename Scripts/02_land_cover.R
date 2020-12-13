library(sf)
library(raster)
#library(rgeos)
#library(ggplot2)
#library(randomForest)
library(tibble)
library(tmap)


# loading data ----
r_sdos <- stack("./raster/LC08_L1TP_224079_20200518_20200527_01_T1_sdos.tif")

posadas <- read_sf('./shp/muni_posadas.shp')

# Processing data ----
# Identifying pxls w/ NA
na_px <- which(!is.na(getValues(r_sdos[[1]])))

# Creating data to Kmeans analysis
mydata <- na.omit(values(r_sdos))
mydata <- scale(mydata) %>% as_tibble() # standardize variables

# training algorithm for unsupervised ----
# K-Means Cluster Analysis
fit <- kmeans(
  mydata, 5, 
  algorithm = c("Hartigan-Wong",
                                "Lloyd",
                                "Forgy",
                                "MacQueen")[1], iter.max = 20)

# Saving statistical results
save(fit, file='./outputs/KmeansFit.RData')

# append cluster assignment
mydata <- data.frame(mydata, fit$cluster) %>% as_tibble()

# creating a raster layer to recieve goup values
landCover <- r_sdos[[1]]
# changing pixel values to group values
landCover[na_px] <- mydata$fit.cluster

map <- tm_shape(landCover) +
  tm_raster(palette = 'cat', style = 'cat', legend.show = FALSE) +
  tm_style(style = 'natural')
tmap_save(map, "./outputs/landcover.png")

writeRaster(landCover, "./raster/landCover.tif", overwrite = TRUE)
