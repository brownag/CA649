library(aqp)
library(soilDB)
library(magrittr)

f <- fetchNASIS()

# estimate PSCS
site(f) <- profileApply(f, function(p) {
    pscs <- estimatePSCS(p)
    data.frame(peiid=profile_id(p), pscs_upr=pscs[1], pscs_lwr=pscs[2])
}, frameify = TRUE)

# calc pscs clay and frags
f.alt <- f %>% 
  glomApply(function(p) c(p$pscs_upr, p$pscs_lwr), truncate = TRUE) %>%
  mutate(thickness = hzdepb - hzdept) %>% 
  mutate_profile(pscs_clay = sum(clay * thickness) / sum(thickness),
                 pscs_frags = sum(fragvoltot * thickness / sum(thickness)))

# transfer pscs calculated values
f$pscs_clay <- f.alt$pscs_clay
f$pscs_frags <- f.alt$pscs_frags

# look at one or more pedons
f.sub <- filter(f, pedon_id %in% c("11PDM013"))

par(mar=c(1,1,1,1))
plotSPC(f.sub, width=0.1, cex.names=0.7, label="pedon_id", id.style="side")

site(f.sub)[,c("pscs_upr","pscs_lwr")]
f.sub$pscs_clay
f.sub$pscs_frags

# depth weighting of quantiles shows heavier textures on average
f.slice <- slice(f, 50:100 ~ clay)
quantile(glomApply(f, function(p) c(50,100), truncate=TRUE)$clay, na.rm=T, probs=c(0.05,0.5,0.95))
quantile(f.slice$clay, na.rm=T, probs=c(0.05,0.5,0.95))

# horizon thickness tends to increase with depth
f$thickness <- f$hzdepb - f$hzdept
f$middepth <- f$hzdept + (f$thickness / 2)
f$thickness[grepl("R|C|BC|qm", f[[hzdesgnname(f)]])] <- NA
f$thickness[is.na(f[[hzdesgnname(f)]])] <- NAha
plot(f$middepth ~ f$thickness, ylim=c(200,0))

hz <- horizons(f)
hz <-  hz[grepl("^A", hz$hzname),]
plot(density(hz$thickness, na.rm=T), main="Horizon Thickness (individual horizons)", xlim=c(0,75))
hz <- horizons(f)
hz <-  hz[grepl("^B", hz$hzname),]
lines(density(hz$thickness, na.rm=T),
      col="RED")
legend("topright", c("A horizons", "B horizons"), col=c("BLACK",'RED'), lty=1)

# calculate cumulative thickness of horizons starting with A
f.a <- f %>% glomApply(function(p) {
  h <- horizons(p)
  h <-  h[grepl("^A", h$hzname),]
  return(c(min(h$hzdept, na.rm=T), max(h$hzdepb, na.rm=T)))
}) %>% mutate_profile(a_hz_thickness = sum(thickness))
plot(density(f.a$a_hz_thickness, na.rm=T), xlim=c(0,75), main="Horizon Thickness (cumulative)")
# calculate cumulative thickness of horizons starting with B
f.b <- f %>% glomApply(function(p) {
  h <- horizons(p)
  h <-  h[grepl("^B", h$hzname),]
  return(c(min(h$hzdept, na.rm=T), max(h$hzdepb, na.rm=T)))
}) %>% mutate_profile(b_hz_thickness = sum(thickness))
lines(density(f.b$b_hz_thickness, na.rm=T),  col="RED")
legend("topright", c("A horizons", "B horizons"), col=c("BLACK",'RED'), lty=1)

horizons(f)[which(f$thickness > 150),]     

filter(f, peiid==1339907)
