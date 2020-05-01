# get project terrain data
library(sf)
library(aqp)
library(soilDB)
library(magrittr)
library(raster)

b <- st_read(".", "project_mapunits_b")

dem <-  raster('C:/Geodata/project_data/MUSum_10m_MLRA/DEM_SON_int_AEA.tif')

b <- sf::st_transform(b, st_crs(dem))
dem.c <- crop(dem, b)
r.b <- fasterize::fasterize(b, dem.c)
dem.b <- mask(dem.c, r.b)
plot(dem.b)

writeRaster(dem.b, filename = "dem_project.tif")
