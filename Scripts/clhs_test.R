library(clhs)
library(raster)
library(rgdal)

# shp <- readOGR(dsn="S:/NRCS/Archive_Andrew_Brown/CA628/20190506_blm.shp")
# shp <- readOGR(dsn="S:/NRCS/Archive_Andrew_Brown/CA649/20190610_HLC.shp")
shp <- readOGR("shp/McClureDonPedro.shp")
# shp <- readOGR("shp/EastBeltMVO.shp")
water <- readOGR('shp/water_I_pbS_1km.shp')
mvo <- readOGR('../../Maps/Mz_metavolcanic_ca.shp')

# if needed: only take biggest feature, after exploding multipart features
#shp <- sp::disaggregate(shp)
#shp <- shp[which(area(shp) == max(area(shp)))[1],]

plot(shp)
# abr <- raster('L:/NRCS/MLRAShared/Geodata/project_data/ssro2_ann_beam_rad_int.tif')
# r <- stack("raster/strahler_6/SAGA_DEM_derived_10m_strahler-level6.tif", "shannon_strahler6.tif")
# abr.sub <- crop(abr, r)
# names(abr.sub) <- 'abr'
# abr.sub <- projectRaster(abr, r)
# r <- stack(r, abr.sub)
# names(r) <- c('elevation','hillshade','slope','aspect','plancurv','profcurv',
#                   'convidx', 'closeddep','catchmentarea', 'twi','ls','channelbase',
#                   'channelvdist', 'valleydepth', 'relativepos',"entropy", "abr")
# writeRaster(r, filename="raster/L6_dem_derived_10m.tif")
r <- stack("raster/L6_dem_derived_10m.tif")
shpt <- spTransform(shp, CRS(proj4string(r)))
watert <- spTransform(water, CRS(proj4string(r)))
mvot <- spTransform(mvo, CRS(proj4string(r)))

r.sub <- crop(r, shpt)
r.sub <- mask(r.sub, shpt)
names(r.sub) <- c('elevation','hillshade','slope','aspect','plancurv','profcurv',
                'convidx', 'closeddep','catchmentarea', 'twi','ls','channelbase',
                'channelvdist', 'valleydepth', 'relativepos',"entropy", "abr")

r.sub$inverse_entropy <- (10 - r.sub$entropy)
r.sub$inverse_entropy <- r.sub$inverse_entropy - min(values(r.sub$inverse_entropy), na.rm=TRUE)
toreplace <- values(r.sub$inverse_entropy)
toreplace[is.na(toreplace)] <- 1E10
values(r.sub$inverse_entropy) <- toreplace

plot(r.sub$inverse_entropy)

# check for layers that are invariant within the study area
#use.idx <- apply(slot(slot(r.sub,'layers'),'values'), 2, sd, na.rm = TRUE) > 0

# use.idx[7] <- FALSE # always remove stream order layer (#7)
# use.idx <- which(!(names(r.sub) %in% c('elevation','hillshade','ls',"closeddep","catchmentarea")))
use.idx <- which((names(r.sub) %in% c('slope','convidx','relativepos','twi','abr','inverse_entropy')))
# use.idx <- as.numeric(which(use.idx))

# only sample the layers with nonzero variance
r.samp <- sampleRegular(r.sub[[names(r.sub)[use.idx]]], 100000, sp = TRUE)
r.samp$owner <- as.character(over(r.samp, shpt)$AGENCY)
r.samp$owner[is.na(r.samp$owner)] <- "non-BLM"
r.samp$is_water <- over(r.samp, watert)$MUSYM == "W"
r.samp$is_water[is.na(r.samp$is_water)] <- FALSE
r.samp$touches_mvo <- as.character(over(r.samp, mvot)$PTYPE)
r.samp$touches_mvo[is.na(r.samp$touches_mvo)] <- "non-MVO"

# constrain sample locations to boundaries of parcel and raster data inputs
r.samp <- r.samp[r.samp$owner == "Bureau of Land Management" & !r.samp$is_water & (r.samp$touches_mvo == "Mzv"),]
r.samp <- r.samp[complete.cases(r.samp@data),]

# select 10 optimized points from the regular grid, conditioned on inverse of information content at grid locations
cl <- clhs(r.samp, 20, cost="inverse_entropy")

plot(cl)
lines(shpt)

writeOGR(cl, dsn = ".", layer = "McClureDonPedro_clhs_points", driver = 'ESRI Shapefile', overwrite_layer = T)
