library(raster)
vdistchannel <- raster('Geodata/derived/vdistchannel.tif')
slope <- tan(raster('Geodata/derived/slope.tif'))*100

lowrelief <- (vdistchannel <= 30)
lowgradient <- (slope <= 30)
lowhills <- (lowrelief + lowgradient)

plot(lowhills)
writeRaster(lowhills == 2, filename="Geodata/derived/lowhills.tif", overwrite=T)

precip <- raster('C:/Geodata/project_data/MUSum_PRISM/final_MAP_mm_800m.tif')
ext <- as(extent(lowhills), 'SpatialPolygons')
crs(ext) <- crs(lowhills)
precip <- crop(precip, spTransform(ext, CRS(proj4string(precip))))

precip.utm <- round(projectRaster(precip, lowhills))
precip.mask <- mask(precip.utm, lowhills)

plot(precip.mask)
writeRaster(precip.mask, "Geodata/derived/precip_mm.tif")

lowprecip1 <- precip.mask <= 18*25.4 # 18" is low cutoff for F200; mid range for 163
lowprecip2 <- precip.mask < 22*25.4  # 22" is low cutoff for F201; high cutoff for 163
lowprecip3 <- precip.mask < 25*25.4  # 25" is high cutoff for F200;

bound <- rgdal::readOGR("Geodata/source/CA649_b.shp")
bound <- spTransform(bound, CRS(proj4string(res1)))

res1 <- lowprecip1 + lowhills == 3
plot(res1, main="Precip < 18\"; Slope < 30%; Vertical dist. from channel < 30m")
lines(bound)

res2 <- (lowprecip2 + lowhills == 3)
plot(res2, main="Precip < 22\"; Slope < 30%; Vertical dist. from channel < 30m")
lines(bound)

res3 <- (lowprecip3 + lowhills == 3)
plot(res3, main="Precip < 25\"; Slope < 30%; Vertical dist. from channel < 30m")
lines(bound)


