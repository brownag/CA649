# get project pedons
library(sf)
library(aqp)
library(soilDB)
library(magrittr)

f <- fetchNASIS()
f.sub <- filter(f, !is.na(f$x_std))
coordinates(f.sub) <- ~ x_std + y_std
proj4string(f.sub) <- "+proj=longlat +datum=WGS84"

ssurgo <- st_read('.', 'project_mapunits_ssurgo')

ssurgo.buf <- ssurgo %>%
  st_union() %>% 
  st_convex_hull()

f.sub.sf <- st_as_sf(as(f.sub, 'SpatialPointsDataFrame')) %>%
              st_transform(st_crs(ssurgo.buf)) %>%
              st_intersection(ssurgo.buf)

plot(ssurgo.buf)
plot(f.sub.sf$geometry, add=T, pch=19, cex=0.5)

table(f.sub.sf$taxonname)

f.mvo <- f.sub %>%
  filter(peiid %in% f.sub.sf$peiid) %>%
  filter(f.sub.sf, grepl("colluvium|residuum", pmkind)) %>%
  filter(f.sub.sf, grepl("metavolcanics|greenstone", pmorigin))

f.mvo.sf <- st_as_sf(as(f.mvo, 'SpatialPointsDataFrame'))

plot(f.mvo.sf$geometry, add=T, pch=19, cex=0.5, col="green")

f.mvo$bedrckdepth <- profileApply(f.mvo, estimateSoilDepth, p="Cr|R|Cd|qm",
                                  no.contact.depth=200, no.contact.assigned=Inf)
f.mvo$no.contact <- !is.finite(f.mvo$bedrckdepth)
plot(density(f.mvo$bedrckdepth, na.rm=T, bw=2.5))

plot(f.mvo %>% filter(no.contact), cex.names=0.75,
     id.style="side", label="pedon_id")

# calculate total surface fragments
f.mvo <- mutate_profile(f.mvo, total_surface_frag = (surface_fgravel + 
                                                     surface_gravel + 
                                                     surface_cobbles + 
                                                     surface_stones + 
                                                     surface_boulders))

# calculate surface horizon
mineral.soil.surface.to.15cm <- function(p) {
  m <- getMineralSoilSurfaceDepth(p)
  return(c(m, m+15))
}

# create site-level attributes from thickest mineral horizon in upper 15cm
f.mvo.surface <- glomApply(f.mvo, mineral.soil.surface.to.15cm, modality="thickest") %>% 
  mutate_profile(surface_texture = texture_class[1], 
                 surface_clay = clay[1])

# add surface clay back into parent SPC for downstream analysis
f.mvo$surface_clay <- f.mvo.surface$surface_clay

# inspect proportional breakdown of texture classes 
table(f.mvo.surface$surface_texture)

# calculate horizons intersecting PSCS
f.mvo.pscs <- glomApply(f.mvo, estimatePSCS) %>% 
  mutate_profile(pscs_clay = weighted.mean(clay, hzdepb - hzdept),
                 pscs_frags = weighted.mean(total_frags_pct, hzdepb - hzdept),
                 clayrat = pscs_clay / surface_clay)

plot(density(f.mvo.surface$surface_clay, na.rm=T, bw=1), col="blue", xlim=c(0,60))
#abline(v=c(18,27,35))
#abline(v=quantile(f.mvo.pscs$surface_clay, na.rm=T), col="blue")
lines(density(f.mvo.pscs$pscs_clay, na.rm=T), col="red")
#abline(v=quantile(f.mvo.pscs$pscs_clay, na.rm=T), col="red")

f.mvo$clayrat <- f.mvo.pscs$clayrat
plot(density(f.mvo.pscs$clayrat, na.rm=T))

ff <- f.mvo %>% 
       filter(clayrat > 3)
f2 <- f.mvo[1:1,]

plotSPC(ff, label="pedon_id", id.style="side", cex.names=1)


plotSPC(f2, label="pedon_id", id.style="side")

rgdal::writeOGR(as(f.mvo, 'SpatialPointsDataFrame'), dsn=".", "MVO_pedons", driver="ESRI Shapefile")


     