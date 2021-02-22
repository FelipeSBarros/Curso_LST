source("./Scripts/ToolBox.R")
rasterPath = "./raster/LC08_L1TP_224079_20201212_20201218_01_T1"
shapePath = './shp/muni_posadas.shp'
demPath = "./raster/dem_30m.tif"

args(pre_processing)
#pre_processing(rasterPath, shapePath, demPath)

# Clasificación con optimization
args(landcover)
landcover(
  rasterPath = "./raster/LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip.tif",
  shapePath = './shp/muni_posadas.shp',
  initialGroups = 3)

automate_suhi(rasterPath = "./raster/LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip.tif",
              ndviLayer = 6,
              ndbiLayer = 10, 
              btLayer = 11, # band temperature layer
              urbanID = 5,
              shapePath = './shp/muni_posadas.shp',
              landcoverPath = './raster/LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip_landCover.tif',
              lstPath = "./raster/LC08_L1TP_224079_20201212_20201218_01_T1/LC08_L1TP_224079_20201212_20201218_01_T1_sdos_clip_lst.tif")

lst_explore()
