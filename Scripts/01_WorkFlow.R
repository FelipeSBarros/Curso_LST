source("./Scripts/tool_box.R")
rasterPath = "./raster/LC08_L1TP_224079_20201212_20201218_01_T1"
shapePath = './shp/muni_posadas.shp'
demPath = "./raster/dem_30m.tif"

args(pre_processing)
pre_processing(rasterPath, shapePath, demPath)
