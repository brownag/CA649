# climate data summarizer by parcel


library(rgdal)
library(raster)

parcel.names <- c("LaPalomaRanch","Hornitos","Juniper","MorningStar","Schoolhouse")

shps <- paste0("S:/NRCS/Archive_Andrew_Brown/CA649/",list("20190610_HLC.shp","20190521_blm.shp","20190522_blm.shp","20190523_blm.shp", "20190524_blm.shp"))
shps <- lapply(shps, readOGR)

cpts <- list("20190610_clhs_points.shp","20190521_clhs_points.shp","20190522_clhs_points.shp","20190523_clhs_points.shp","20190524_clhs_points.shp")
cpts <- lapply(as.list(1:5), function(x) {
  s <- readOGR(cpts[[x]])
  s$ident <- paste0(parcel.names[x],"_",0:(length(s)-1))
  return(s)
})
cpts.all <- do.call('rbind', lapply(cpts, function(y) { 
  y[,c(1,length(names(y)))] 
  }))
names(cpts.all) <- c("slope","ident")
cpts.all$slope <- tan(cpts.all[[1]])*100
writeOGR(cpts.all, dsn=".",layer="clhs_sets_1-5", driver="ESRI Shapefile", overwrite_layer = T)

demd <- stack("ca649_terrain_derivatives_10m.tif")

covariate.names <- c("slope","aspect","cross_curv","long_curv","conv_idx",
                     "closed_depression","flow_acc","TWI","LS","channel_base",
                     "vdist_channel","valley_depth","relative_slope_pos",
                     "TRI","MRVBF","MRRTF","VRM")
names(demd) <- covariate.names

prism <- list("L:/NRCS/MLRAShared/Geodata/project_data/MUSum_PRISM/final_MAAT_800m.tif",
              "L:/NRCS/MLRAShared/Geodata/project_data/MUSum_PRISM/final_MAP_mm_800m.tif",
              "L:/NRCS/MLRAShared/Geodata/project_data/MUSum_PRISM/effective_precipitation_800m.tif",
              "L:/NRCS/MLRAShared/Geodata/project_data/MUSum_PRISM/ffd_50_pct_800m.tif",
              "L:/NRCS/MLRAShared/Geodata/project_data/MUSum_PRISM/ffd_90_pct_800m.tif")
prism <- stack(prism)

grds <- lapply(shps, function(s) {
  s2 <- spTransform(s, CRS(proj4string(demd)))
  demc <- crop(demd, s2)
  return(sampleRegular(demc, s$acres*2))
})

clim <- do.call('rbind',lapply(shps, function(s) {
  s2 <- spTransform(s, CRS(proj4string(prism)))
  clic <- crop(prism, s2)
  return(colMeans(sampleRegular(clic, s$acres / 10)))
}))
row.names(clim) <- parcel.names
knitr::kable(clim[order(clim[,1], decreasing = T),])

idx <- 1
lapply(grds, function(g) {
  plot(density(tan(g[,1]) * 100, na.rm=T), main=paste0(parcel.names[idx], " Slope"))
  clhs.slopes <- tan(cpts[[idx]][[1]]) * 100
  points(clhs.slopes, rep(0, length(clhs.slopes)))
  idx <<- idx + 1 # >:)
})


