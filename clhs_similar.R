# clhs - similar location finder
library(clhs)
library(raster)
library(rgdal)

setwd("E:/CA649")
# shp <- readOGR(dsn="S:/NRCS/Archive_Andrew_Brown/CA628/20190506_blm.shp")
shp <- readOGR(dsn="S:/NRCS/Archive_Andrew_Brown/CA649/20190610_HLC.shp")

# only take biggest feature, after exploding multipart features
shp <- sp::disaggregate(shp)
shp <- shp[which(area(shp) == max(area(shp)))[1],]

r <- stack("E:/CA649/ca649_terrain_derivatives_10m.tif")
shpt <- spTransform(shp, CRS(proj4string(r)))

r.sub <- crop(r, shpt)
r.sub <- mask(r.sub, shpt)

# check for layers that are invariant within the study area
use.idx <- apply(slot(slot(r.sub,'data'),'values'), 2, sd, na.rm = TRUE) > 0

use.idx[7] <- FALSE # always remove stream order layer (#7)
use.idx <- as.numeric(which(use.idx))

# only sample the layers with nonzero variance
r.samp <- sampleRegular(r.sub[[names(r.sub)[use.idx]]], 100000, sp = TRUE)
clhs.all <- readOGR(dsn = ".", layer = "20190610_clhs_points")

library(Hmisc)
exhaustive <- r.samp[complete.cases(r.samp@data),]
exhaustive$source <- "exhaustive"

clhs.set <- clhs.all[,1:16]
names(clhs.set) <- names(r.samp)
clhs.set$source <- "cLHS"

dat <- rbind(exhaustive, clhs.set)
#dmat <- dist(dat@data[,1:16])
#hc <- hclust(dmat)
#memb <- cutree(hc, k=5)

hc <- cluster::clara(dat@data[,1:16], 5)
memb <- hc$clustering
dat$member <- memb
rdat <- rasterize(exhaustive, raster(extent(exhaustive), res=10), dat$member, fun=min)
proj4string(rdat) <- proj4string(exhaustive)

plot(rdat, col=viridis::viridis(5))
points(clhs.all, pch=24, cex=1, col="yellow")
lines(shpt, col="WHITE")

writeRaster(rdat, "LaPalomaHRC_cLHS_membership.tif", overwrite=T)
