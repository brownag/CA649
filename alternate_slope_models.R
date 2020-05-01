library(raster)
r <- raster("E:/CA649/dem_project.tif")
r.c <- crop(r, extent(r))
s_queen <- terrain(r.c, unit="radians")
s_rook <- terrain(r.c, unit="radians", neighbors=4)
plot(density(abs(values(s_queen) - values(s_rook)), na.rm=T))

writeRaster(s_queen, filename="E:/CA649/slope_queen_project.tif", overwrite=T)
writeRaster(s_rook, filename="E:/CA649/slope_rook_project.tif", overwrite=T)

s_queen.c <- crop(s_queen, extent(s_queen) / 100)

par(mfrow=c(1,2))
plot(s_queen.c)
plot(focal.max(s_queen.c,3))

focal.max <- function(r, n, ...) {
  focal(r, w=matrix(1,n,n), fun=max, ...)
}

writeRaster(focal.max(s_queen, 3) / (2*pi) * 100, 
            filename="E:/CA649/slopemax_queen_project.tif", 
            overwrite=T)

writeRaster(focal.max(s_rook, 3) / (2*pi) * 100, 
            filename="E:/CA649/slopemax_rook_project.tif", 
            overwrite=T)

focal.q <- function(r, n, q, ...) {
  focal(r, w = matrix(1, n[1], n[1]), 
        fun = function(x) {
          as.numeric(quantile(x, probs=q[1], na.rm=TRUE))
        }, ...)
}

par(mfrow=c(1,3))
plot(s_queen.c)
plot(focal.max(s_queen.c, 3))
plot(focal.q(s_queen.c, 3, 0.5))

res <- focal.q(s_queen, 3, 0.8, filename="foo.tif")
writeRaster(res * 2*pi, 
            filename="E:/CA649/slope_q80_queen_project.tif", 
            overwrite=T)
