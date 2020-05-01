library(aqp)
library(soilDB)
library(raster)

f <- fetchNASIS()
k <- fetchKSSL(mlra="18")

prism_precip <- raster("C:/Geodata/project_data/MUSum_PRISM/final_MAP_mm_800m.tif")

coordinates(k) <- ~ x + y
proj4string(k) <- "+proj=longlat +datum=WGS84"

k$prism_map_mm <- extract(prism_precip, slot(k, 'sp'))

# all MLRA 18 ranges
lhz <- horizons(k)
quantile(lhz[grepl("^A", lhz$hzn_desgn),]$estimated_om, na.rm=TRUE)

k.sub1 <- filter(k, prism_map_mm < 21*25.4)
omgt5 <- filter(k.sub1, !(estimated_om > 20))
om_slab <- slab(slice(filter(omgt5, !pedon_key %in% c("20200")), 0:25 ~ estimated_om), fm = ~estimated_om)
om_slab$mid <- (om_slab$top + om_slab$bottom) / 2
plot(om_slab$mid ~ om_slab[,c(5)], type="l", ylim=c(25,0), xlim=c(0,10), lwd=2)
lines(om_slab$mid ~ om_slab[,c(4)], lty=2, lwd=2)
lines(om_slab$mid ~ om_slab[,c(6)], lty=2, lwd=2)
lines(om_slab$mid ~ om_slab[,c(3)], lty=2, lwd=1)
lines(om_slab$mid ~ om_slab[,c(7)], lty=2, lwd=1)
abline(v=c(1.5, 3.5, 5.5))
abline(v=c(1, 2, 3), col="green")

omgt5$theme <- factor(omgt5$estimated_om > 5)
omgt5$taxonname <- factor(omgt5$taxonname)
omgt5$estimated_om <- round(omgt5$estimated_om, 2)

groupedProfilePlot(omgt5, label="pedon_id", name="estimated_om",
                   color="theme", groups="taxonname")

lhz1 <- horizons(k.sub1)
quantile(lhz1[grepl("^A", lhz1$hzn_desgn),]$estimated_om, probs=c(0.05,0.5,0.95),na.rm=TRUE)

plot(density(lhz1[grepl("^A", lhz1$hzn_desgn),]$estimated_om, na.rm=TRUE))
lines(density(lhz[grepl("^A", lhz$hzn_desgn),]$estimated_om, na.rm=TRUE), col='blue')
