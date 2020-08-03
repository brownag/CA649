# soil aspect model (using beam radiance)
#  this investigation follows from the determination that abr is an important predictor in soil depth

aspect.phases <- split.points(f, 'abr', c(0,60000,71000), c(60000,71000,100000))
f$cool.phase <- as.logical(aspect.phases[,1])
f$unphased <- as.logical(aspect.phases[,2])
f$warm.phase <- as.logical(aspect.phases[,3])

groupedProfilePlot(filter(f, cool.phase), groups="taxonname")

groupedProfilePlot(filter(f, unphased, taxonname == "Dunstone"), groups="taxonname")
abline(h=50, lwd=2, col="red")

groupedProfilePlot(filter(f, warm.phase), groups="taxonname")

warm <-  filter(f, warm.phase)
cool <- filter(f, cool.phase)

plot(density(f$slope_field, na.rm=T, bw=5), ylim=c(0,0.04))
lines(density(warm$slope_field, bw=5), col="red")
lines(density(cool$slope_field, bw=5, na.rm=T), col="blue")

plot(density(f$bedrckdepth, na.rm=T, bw=5), ylim=c(0,0.04))
lines(density(warm$bedrckdepth, bw=5), col="red")
lines(density(cool$bedrckdepth, bw=5, na.rm=T), col="blue")
abline(v=median(f$bedrckdepth, na.rm=T), col="black")
abline(v=median(warm$bedrckdepth), col="red")
abline(v=median(cool$bedrckdepth), col="blue")

plot(data=site(cool), bedrckdepth ~ slope_field)

plot(slot(filter(f, warm.phase),'sp'), col='red')
points(slot(filter(f, cool.phase),'sp'), col='blue', pch="+")
points(slot(filter(f, !warm.phase, !cool.phase),'sp'), pch=".", cex=1.5)

b <- rgdal::readOGR("Geodata/source/CA649_b.shp")
b <- sp::spTransform(b, CRS(proj4string(f)))
lines(b)

pb <- rgdal::readOGR("project_mapunits_b.shp")
b <- rgdal::readOGR("Geodata/source/CA649_b.shp")
plot()
rgdal::writeOGR(crop(pb, extent(b)),dsn=".", "mariposa_project_b", driver="ESRI Shapefile")
