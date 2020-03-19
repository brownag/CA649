# get spatial data from project mapunit list
library(rgdal)
library(soilDB)

project.mapunits <- read.csv('project_mapunits.csv')
s <- fetchSDA_spatial(project.mapunits$nationalmusym,
                 by.col="nmusym",
                 add.fields=c('mapunit.musym','mapunit.muname'),
                 chunk.size = 1)
plot(s)
s <- spTransform(s, CRS("+proj=utm +zone=11 +datum=WGS84"))
writeOGR(s, ".", "project_mapunits_ssurgo", driver="ESRI Shapefile", overwrite_layer = TRUE)
