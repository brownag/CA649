# polygon QC script for vertex density
# identify polygons/lines with abnormally high or low vertex density

library(sf)
library(dplyr)

# load "working progress" file -- combined project mapunit extent CA630+CA649
wp <- st_read("working_progress.shp")

# take just the 649 delineations
#wp <- wp[which(is.na(wp$arsymbl) | wp$arsymbl != "CA630"),]

# calculate vertex density (meters per vertex [by polygon])
wp$poly.vert.density <-  lwgeom::st_perimeter(wp) / mapview::npts(st_geometry(wp), by_feature = TRUE)

# inspect
plot(density(wp$poly.vert.density))

# determine quantiles
quantile(wp$poly.vert.density, probs=c(0,0.01,0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.99,1))

# apply quantile-based threshold for output (>90th percentile or so)
check.polys <- wp[which(as.numeric(wp$poly.vert.density) > 35),]

# spatial plot
plot(check.polys$geometry, col=check.polys$poly.vert.density)

# more detailed: now, lets evaluate subsets to identify bad portions of lines
# convert to LINESTRING
wpl <- st_cast(st_cast(wp, 'MULTILINESTRING'), 'LINESTRING')

# convert LINESTRING to constituent POINTs
linepoints <- lapply(wpl$geometry, function(x) {
  st_sfc(x) %>% st_cast(., "POINT")
})

# split each poly feature (single LINESTRING) into _multiple_ LINESTRING
#   use an arbitrary number of vertices -- 30 should be enough to calculate local density
linesplit <- sapply(linepoints, function(p) {
  # make n chunks of size 30 for each polygon
  idx <- soilDB::makeChunks(1:length(p), size=30)
  
  buf <- list()
  for(i in 1:max(idx)) {
    # subset vertices by chunk
    p.sub <- p[which(idx == i)]
    
    # short circuit
    if(!length(p.sub))
      return(NULL)
    
    # inspect
    # plot(st_cast(st_combine(p.sub), 'LINESTRING'))
    
    # recombine subset of POINT (vertices) into LINESTRING
    buf[[i]] <- st_cast(st_combine(p.sub), 'LINESTRING')
  }
  return(buf)
})

# unpack results and combine into sf object
lines_qc <- st_sf('geometry'=do.call('c', do.call('c', linesplit)), crs=st_crs(wp))

# meters per vertex (by line)
lines_qc$line.vert.density <-  st_length(lines_qc) / mapview::npts(st_geometry(lines_qc), by_feature = TRUE)

plot(density(lines_qc$line.vert.density, from=0, bw = 5))
quantile(lines_qc$line.vert.density, probs=c(0,0.01,0.05,0.10,0.25,0.5,0.75,0.9,0.95,0.99,1))

# apply threshold for output (biggest 1-2% vertex spacings)
check.lines <- lines_qc %>% 
  filter(line.vert.density > units::as_units(50, "m"))

plot(check.lines, col=round(check.lines$line.vert.density))

# polygons with one or more "bad" vertices
idx <- which(apply(as.matrix(st_intersects(wp, check.lines)), 1, any))
polys.badlines <- wp[idx,]
plot(polys.badlines$geometry, col="red")

st_write(check.polys, "bad_polygons.shp", delete_layer = TRUE)

st_write(check.lines, "bad_lines.shp", delete_layer = TRUE)

st_write(polys.badlines, "bad_lines_poly.shp", delete_layer = TRUE)
