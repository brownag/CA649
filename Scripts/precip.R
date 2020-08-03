
precip <- raster('C:/Geodata/project_data/MUSum_PRISM/final_MAP_mm_800m.tif')
ext <- as(extent(lowhills), 'SpatialPolygons')
crs(ext) <- crs(lowhills)
precip <- crop(precip, spTransform(ext, CRS(proj4string(precip))))

precip.utm <- round(projectRaster(precip, lowhills))
precip.mask <- mask(precip.utm, lowhills)

plot(precip.mask)
writeRaster(precip.mask, "Geodata/derived/precip_mm.tif")