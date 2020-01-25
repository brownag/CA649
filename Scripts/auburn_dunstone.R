library(aqp)
library(soilDB)

# this script determines whether a soil has paralithic, paralithic over lithic, or lithic
# or indeterminate (incompletely described). 

# by default, paralithic materials _exactly_ 25cm thick are indeterminate.

series <- c("Auburn","Bonanza","Dunstone","Exchequer","Millerton")

# load all CA630 pedons
f <- fetchNASIS()
f <- f[!is.na(f$x_std),]
coordinates(f) <- ~ x_std + y_std
proj4string(f) <- "+proj=longlat +datum=WGS84"


f.sub <- f[f$taxonname %in% series,]
f.geom <- SDA_spatialQuery(as(f.sub, 'SpatialPointsDataFrame'), what="geom")

lithic.depth <- profileApply(f.sub, function(p) {
  p$hzdept[grep(p$hzname, pattern="R", ignore.case = FALSE)[1]]
})

paralithic.depth <- profileApply(f.sub, function(p) {
  p$hzdept[grep(p$hzname, pattern="Cr", ignore.case = FALSE)[1]]
})

max.description.depth <- profileApply(f.sub, function(p) {
  max(p$hzdepb)
})

lithic.depth[(lithic.depth == max.description.depth)] <- NA
paralithic.depth[(paralithic.depth == max.description.depth)] <- NA

f.sub$cr_thickness <- profileApply(f.sub, function(p) {
  idx <- grep(p$hzname, pattern="Cr")
  return(sum(p$hzdepb[idx] - p$hzdept[idx]))
})

f.sub$indeterminate <- (is.na(lithic.depth) & f.sub$cr_thickness <= 25)
f.sub$is_lithic <- is.na(lithic.depth - paralithic.depth) & !is.na(lithic.depth) & !f.sub$indeterminate
f.sub$is_paralithic <- is.na(lithic.depth - paralithic.depth) & is.na(lithic.depth) & !f.sub$indeterminate
f.sub$is_paraoverlithic <- !is.na(lithic.depth - paralithic.depth) & ((lithic.depth - paralithic.depth) <= 25) & !is.na(lithic.depth) & !f.sub$indeterminate

par(mar=c(0,0,0,0))
plot(aggregate(f.geom))
points(as(f.sub[f.sub$indeterminate],'SpatialPointsDataFrame'), pch=23, cex=1.5, bg="BLACK")
points(as(f.sub[f.sub$is_lithic],'SpatialPointsDataFrame'), pch=23, cex=1.5, bg="RED")
points(as(f.sub[f.sub$is_paralithic],'SpatialPointsDataFrame'), pch=23, cex=1.5, bg="BLUE")
points(as(f.sub[f.sub$is_paraoverlithic],'SpatialPointsDataFrame'), pch=23, cex=1.5, bg="GREEN")
legend("topright",legend=paste0(c("Indeterminate","Lithic", "Paralithic", "Paralithic over Lithic"), " (n = ", 
                                  colSums(site(f.sub)[,c("indeterminate","is_lithic","is_paralithic","is_paraoverlithic")]), ")"), 
       pch=NA,
       cex=2,
       fill=c("BLACK","RED","BLUE","GREEN"))

plot(f.sub[f.sub$is_lithic])
plot(f.sub[f.sub$is_paralithic])
plot(f.sub[f.sub$is_paraoverlithic])
plotSPC(f.sub[f.sub$indeterminate], max.depth=200, label = 'pedon_id')
