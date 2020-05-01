# get project pedons
library(sf)
library(aqp)
library(soilDB)
library(magrittr)
library(raster)

f <- fetchNASIS()
f.sub <- filter(f, !is.na(x_std))

coordinates(f.sub) <- ~ x_std + y_std
proj4string(f.sub) <- "+proj=longlat +datum=WGS84"

ssurgo.buf <- st_read('.', "project_mapunits_b")

f.sub.sf <- st_as_sf(as(f.sub, 'SpatialPointsDataFrame')) %>%
  st_transform(st_crs(ssurgo.buf)) %>%
  st_intersection(ssurgo.buf)

prism.map <- raster('C:/Geodata/project_data/MUSum_PRISM/final_MAP_mm_800m.tif')
prism.maat <- raster('C:/Geodata/project_data/MUSum_PRISM/final_MAAT_800m.tif')
terrain.slope <- raster("C:/Geodata/project_data/MUSum_10m_MLRA/Slope_SON_int_AEA.tif")

plot(ssurgo.buf$geometry)
plot(f.sub.sf$geometry, add=T, pch=19, cex=0.5)

table(f.sub.sf$taxonname)

f.mvo <- f.sub %>%
  filter(peiid %in% f.sub.sf$peiid)# %>%
  #filter(f.sub.sf, grepl("colluvium|residuum", pmkind)) %>%
  #filter(f.sub.sf, grepl("metavolcanics|greenstone", pmorigin))

f.mvo$map <- raster::extract(prism.map, f.mvo.sf)
f.mvo$maat <- raster::extract(prism.maat, f.mvo.sf)

# anything with a C in it is considered substratum
# several pedons are in weathered bedrock ~thick paralithic VD soils, but not that common
f.mvo$bedrckdepth <- profileApply(f.mvo, estimateSoilDepth, p="^[23B]*C|R|Cd|qm",
                                  no.contact.depth=150, no.contact.assigned=Inf)
f.mvo$no.contact <- !is.finite(f.mvo$bedrckdepth)
plot(density(f.mvo$bedrckdepth, na.rm=T, bw=2.5))

# this one is a big ol red and deep one
plot(f.mvo %>% filter(no.contact), cex.names=0.75,
     id.style="side", label="pedon_id")

# calculate total surface fragments
f.mvo <- mutate_profile(f.mvo, total_surface_frag = (surface_fgravel + 
                                                       surface_gravel + 
                                                       surface_cobbles + 
                                                       surface_stones + 
                                                       surface_boulders))
# mvo peiids
peiid.mvo <- paste0(as.numeric(f.mvo$peiid),collapse=",")

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
f.mvo.pscs <- glomApply(f.mvo, estimatePSCS, truncate=TRUE) %>% 
  mutate_profile(pscs_clay = weighted.mean(clay, hzdepb - hzdept),
                 pscs_frags = weighted.mean(total_frags_pct, hzdepb - hzdept),
                 clayrat = pscs_clay / surface_clay,
                 claymax = max(clay, na.rm=T))

plot(density(f.mvo.surface$surface_clay, na.rm=T, bw=1), col="blue", xlim=c(0,60))
#abline(v=c(18,27,35))
#abline(v=quantile(f.mvo.pscs$surface_clay, na.rm=T), col="blue")
lines(density(f.mvo.pscs$pscs_clay, na.rm=T), col="green")
lines(density(f.mvo.pscs$claymax, na.rm=T), col="red")
#abline(v=quantile(f.mvo.pscs$pscs_clay, na.rm=T), col="red")

f.mvo$clayrat <- f.mvo.pscs$clayrat
plot(density(f.mvo.pscs$clayrat, na.rm=T))

# look at soils with high 'eluviation ratio'
ff <- f.mvo %>% 
  filter(clayrat > 2.25)
plotSPC(ff, label="pedon_id", id.style="side", cex.names=1)

# calculate pscs_clay, pscs_frags, clay ratio and maximum clay
f.mvo.pscs <- glomApply(f.mvo, estimatePSCS, truncate = TRUE) %>% 
  mutate_profile(pscs_clay = weighted.mean(clay, hzdepb - hzdept),
                 pscs_frags = weighted.mean(total_frags_pct, hzdepb - hzdept),
                 clayrat = pscs_clay / surface_clay,
                 claymax = max(clay, na.rm=T))

# merge back into site table
site(f.mvo) <- site(f.mvo.pscs)[,c(idname(f.mvo),"pscs_clay","pscs_frags")]

# write site table to shapefile
rgdal::writeOGR(as(f.mvo, 'SpatialPointsDataFrame'), 
                dsn=".", "MVO_pedons", 
                driver="ESRI Shapefile", overwrite_layer = TRUE)

very.shallows <- f.mvo.pscs %>% 
  filter(psctopdepth == 0)

matchSPC <- function(spc_in, spc_out) {
  i1 <- idname(spc_out)
  i2 <- idname(spc_in)
  
  if(i1 != i2)
    warning("match: profile ID names do not match (",i1,":",i2,")", call. = FALSE)
  
  idx <- which(site(spc_out)[[idname(spc_out)]] %in% spc_in[[idname(spc_in)]])
  
  if(!length(idx))
    idx <- 0
  
  return(spc_out[idx,])
}

deeper <- f.mvo.pscs %>% 
  filter(psctopdepth > 0)
plot(deeper)


f.mvo.sub <- f.mvo %>% 
  glomApply(.fun=function(x) return(c(0,15)), modality="thickest") %>%
  filter(texcl == "cl") %>%
  match(f.mvo)

plot(f.mvo.sub)

hz.thickness <- function(p, hzdesgn) {
  h <- horizons(p)
  d <- horizonDepths(p)
  idx <- grep(hzdesgn, h[[hzdesgnname(p)]])
  if(length(idx)) {
    thk <- h[idx, d[2]] - h[idx, d[1]]
    return(sum(thk))
  } else { return(0) }
}

f.mvo$cr_thickness <- profileApply(f.mvo, hz.thickness, hzdesgn="Cr")
f.mvo$r_thickness <- profileApply(f.mvo, hz.thickness, hzdesgn="R")
f.mvo$bt_thickness <- profileApply(f.mvo, hz.thickness, hzdesgn="t")

# factor taxonname for plotting and summaries
f.mvo$taxonname <- factor(f.mvo$taxonname)

## get just 163/200 soils
f.mvo.low <- f.mvo %>%
  filter(map <= 21*25.4)

## get higher elevation band
f.mvo.hi <- f.mvo %>%
  filter(map > 21*25.4)

# taxa by precip band
groupedProfilePlot(f.mvo.low, groups="taxonname")
res <- table(f.mvo.low$taxonname)
res <- res[rev(order(res))]
res

# fine psc
groupedProfilePlot(filter(f.mvo.low, pscs_frags >= 35), groups="taxonname")

groupedProfilePlot(f.mvo.hi, groups="taxonname")     
res <- table(f.mvo.hi$taxonname)
res <- res[rev(order(res))]
res

# inspect some soils that have major diagnostic features identified
groupedProfilePlot(f.mvo %>% 
                     filter(cr_thickness < 10, cr_thickness > 0, r_thickness > 0, bt_thickness > 0), 
                   label = "pedon_id", id.style = "side", cex.names=0.8, groups="taxonname")

# let's just set up NASIS for the low elevation soils
plot(f.mvo.hi@sp)
paste0(f.mvo.hi$peiid, collapse=",")

f.all <- f.mvo.hi
