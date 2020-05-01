# get spatial data from project mapunit list
library(sf)
library(rgdal)
library(soilDB)

p.names <- list("MLRA 18 - CA649 Auburn & Blasingame 999 Map Units on Volcanic Rocks",
                "MLRA 18 - Trabuco TAD2 & TaE2 Mariposa Area (CA649)",
                "MLRA 18 - CA649 Metavolcanic Rockland (PLACEHOLDER)",
                "MLRA 18 - Metavolcanic MLRA Mapunits (PLACEHOLDER)",
                "MLRA 18 - CA649 Metavolcanic pESD Mapunits (PLACEHOLDER)")

project.mapunits <- do.call('rbind', lapply(p.names, get_projectmapunit_from_NASISWebReport))

write.csv(project.mapunits, file = 'project_mapunits.csv')

# s <- fetchSDA_spatial(res$nationalmusym,
#                  by.col="nmusym",
#                  add.fields=c('mapunit.musym','mapunit.muname'),
#                  chunk.size = 1)
# ssurgo <- spTransform(s, CRS("+proj=utm +zone=11 +datum=WGS84"))
# rgdal::writeOGR(ssurgo, ".", "project_mapunits_ssurgo", driver="ESRI Shapefile", overwrite_layer = TRUE)

ssurgo <- rgdal::readOGR(".", "project_mapunits_ssurgo")

ssurgo.buf <- st_as_sf(ssurgo) %>%
  st_union() %>% 
  st_convex_hull()

st_write(ssurgo.buf, dsn = "project_mapunits_b.shp", update=TRUE)

newgdb <- st_read('Geodata/Offical_Geodatabase/FGCA649_Projects_2020_0323_agb.gdb', 'ca649_a')

ssurgo <- st_read('.', 'project_mapunits_ssurgo')
ssurgo <- st_transform(ssurgo, st_crs(newgdb))

ca649.mukeys <- project.mapunits$lmapunitiid[project.mapunits$areasymbol == "CA649"]
ssurgo.no649 <- ssurgo[!ssurgo$mukey %in% ca649.mukeys, 'musym']
ssurgo.no649$orig_musym <- ssurgo.no649$musym
ssurgo.no649 <- merge(ssurgo.no649, project.mapunits[project.mapunits$areasymbol == "CA630",])
working.musyms <- c(project.mapunits$musym,
                    '7083b','7085b','7076b','7078b','7079b','8033','8034')

# prelim low elevation MUs
newgdb.pmu <- newgdb[newgdb$MUSYM %in% working.musyms, c('MUSYM', 'orig_musym')]
names(newgdb.pmu) <- c("musym", "orig_musym", "geometry")
st_geometry(newgdb.pmu) <- "geometry"
#newgdb.pmu <- merge(newgdb.pmu, project.mapunits, all.x=TRUE)

res <- rbind(ssurgo.no649[,names(newgdb.pmu)[1:2]], newgdb.pmu)
plot(res$geometry)

res$mlrassoarea <- NULL
res$nonmlrassaarea <- NULL
res$projecttypename <- NULL
res$fiscalyear <- NULL
res$areaname <- NULL
res$pmu_seqnum <- NULL
head(res)

res2 <- merge(res, project.mapunits[project.mapunits$areasymbol == "CA630",], all.x=TRUE)
plot(res2$geometry)

res3 <- st_cast(res2, "POLYGON")

st_write(res3, "working_progress.shp", delete_layer = TRUE, update=TRUE)
library(dplyr)
update_impact <- res3 %>% 
  dplyr::filter(musym %in% c("7083b","7085b","7076b","7078b","7079b",
                             "7083","7085","7076","7078","7079","8033","8034") &
                is.na(areasymbol))

plot(update_impact$geometry, col="BLUE")

update_impact_b <- st_union(update_impact)
plot(update_impact_b)

original_ssurgo <- st_read("Geodata/Offical_Geodatabase/FGCA649_Projects_2019_1125_ra.gdb", layer = "ca649_a")
ssurgo_updated <- st_intersection(original_ssurgo, update_impact_b)
plot(ssurgo_updated$Shape)

res <- ssurgo_updated %>% 
  group_by(MUSYM) %>% 
  dplyr::summarize(acreage=units::set_units(sum(st_area(Shape)), acre))
View(res)
