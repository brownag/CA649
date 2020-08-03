# clhs bootstrap density

# run clhs routine exhaustively, and analyze density of optimal point locations
n.target <- 10
n.iterations <- 10000

setwd("E:/CA649")
# shp <- readOGR(dsn="S:/NRCS/Archive_Andrew_Brown/CA628/20190506_blm.shp")
shp <- readOGR(dsn="S:/NRCS/Archive_Andrew_Brown/CA649/20190610_HLC.shp")

# only take biggest feature, after exploding multipart features
shp <- sp::disaggregate(shp)
shp <- shp[which(area(shp) == max(area(shp)))[1],]

r <- stack("E:/CA649/ca649_terrain_derivatives_10m.tif")
shpt <- spTransform(shp, CRS(proj4string(r)))
r.sub <- crop(r, shpt)
r.sub.mask <- raster::mask(r.sub, shpt)

# check for layers that are invariant within the study area
use.idx <- apply(slot(slot(r.sub.mask, 'data'), 'values'), 2, sd, na.rm = TRUE) > 0
use.idx[7] <- FALSE # always remove stream order layer (#7)
use.idx <- as.numeric(which(use.idx))

# only sample the layers with nonzero variance
r.samp <- sampleRegular(r.sub.mask[[names(r.sub.mask)[use.idx]]], 
                        size = 1000, sp = TRUE, na.rm = TRUE)
r.samp <- r.samp[which(complete.cases(over(r.samp, shpt))),]

# select n.target optimized points n.iterations times (with 100 iterations per sample set)
res <- lapply(as.list(1:n.iterations), function(x) { clhs(r.samp, n.target, iter=2, progress=F) })
res.all <- do.call('rbind', res)
res.all$ID <- 1

library(ggmap)
register_google(key="AIzaSyCCfQVrzF93WcS59cegFTUiZxIB7fd9350")

res.all.t <- spTransform(res.all[complete.cases(res.all.t@data),], CRS("+proj=longlat +datum=WGS84"))
my.coords <- coordinates(res.all.t)
colnames(my.coords) <- c("lon","lat")
mydaty <- cbind(my.coords, res.all.t@data)

basemap <- get_map(location = c(mean(extent(res.all.t)[1:2]),
                                mean(extent(res.all.t)[3:4])), zoom = 14)

clhs.actual <- readOGR(dsn = ".", layer = "20190610_clhs_points")
clhs.actual <- spTransform(clhs.actual, CRS("+proj=longlat +datum=WGS84"))
clhs.locs <- coordinates(clhs.actual)
colnames(clhs.locs) <- c("lon","lat")
clhs.points <- cbind(clhs.locs, clhs.all@data)

plot(raster(MASS::kde2d(my.coords[,1], my.coords[,2])))

ggmap(basemap) + 
  stat_density2d(aes(fill = ..level..), alpha = .25, geom = "polygon", data = mydaty) + 
  scale_fill_viridis() +
  geom_point(data=clhs.points)

