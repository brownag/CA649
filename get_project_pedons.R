# get project pedons
library(sf)
library(aqp)
library(soilDB)
library(magrittr)

f <- fetchNASIS()
f.sub <- filter(f, !is.na(f$x_std))
coordinates(f.sub) <- ~ x_std + y_std
proj4string(f.sub) <- "+proj=longlat +datum=WGS84"

ssurgo <- st_read('.', 'project_mapunits_ssurgo')

f.sub.sf <- st_as_sf(as(f.sub, 'SpatialPointsDataFrame')) %>%
              st_transform(st_crs(ssurgo)) %>%
              st_intersection(ssurgo)
plot(f.sub.sf)
