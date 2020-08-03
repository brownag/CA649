library(aqp)
library(soilDB)
library(sf)

# read individual metavolcanic formation boundaries
jgo <- st_read("Geodata/derived/metavolcanics_Jgo_south.shp")
jpb <- st_read("Geodata/derived/metavolcanics_Jpb.shp")
jmvo <- rbind(jgo, jpb)

# get selected set
f <- fetchNASIS()
f <- filter(f, !is.na(x_std), x_std < -115)
coordinates(f) <- ~ x_std + y_std
proj4string(f) <- "+proj=longlat +datum=WGS84"
plot(slot(f,'sp'))

# calculate which mvo feature each point is in
pit.locations <- st_as_sf(slot(f, 'sp'))
pit.locations$profile_id <- profile_id(f)
pit.locations <- st_transform(pit.locations, st_crs(jmvo))
pit.locations$poly <- unlist(lapply(st_intersects(pit.locations, jmvo), function(x) {
  if(!length(x)) {
    return(0) 
  }
  x
}))
f$poly <- pit.locations$poly
plot(slot(f, 'sp'), col=f$poly)

jpb.pedons <- filter(f, poly %in% c(4,5))
