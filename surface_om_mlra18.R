library(aqp)
library(soilDB)

mlra18 <- fetchKSSL(mlra="18")
coordinates(mlra18) <- ~ x + y
proj4string(mlra18) <- "+proj=longlat +datum=WGS84"

# use glom to get upper 10cm of mineral soil from each profile
upper10mineral <- function(spc) { 
  glomApply(spc, function(p) {
    top <- getMineralSoilSurfaceDepth(p, hzdesgn=hzdesgnname(p))
    bottom <- top + 10
    return(c(top, bottom))
  }, truncate=TRUE)
}
mlra18.upper10mineral <- upper10mineral(mlra18)

# combine some classes -- to simplify
lcos.idx <- grep("LCOS", mlra18.upper10mineral$lab_texture_class, ignore.case = TRUE)
sandy.loam.idx <- grep("SL", mlra18.upper10mineral$lab_texture_class, ignore.case = TRUE)
sic.idx <- grep("SIC", mlra18.upper10mineral$lab_texture_class, ignore.case = TRUE)

mlra18.upper10mineral$lab_texture_class_simp <- toupper(mlra18.upper10mineral$lab_texture_class)
mlra18.upper10mineral$lab_texture_class_simp[lcos.idx] <- "LS" 
mlra18.upper10mineral$lab_texture_class_simp[sandy.loam.idx] <- "SL" 
mlra18.upper10mineral$lab_texture_class_simp[sic.idx] <- NA

plot(data=horizons(mlra18.upper10mineral), estimated_om ~ factor(lab_texture_class_simp),
     main='Organic Matter % - upper 10cm - MLRA 18', ylab="Organic Matter %",
     xlab="Fine-earth Fraction Texture")

plot(data=horizons(mlra18.upper10mineral), estimated_om ~ clay,
     main='Organic Matter % - upper 10cm - MLRA 18', ylab="Organic Matter %",
     xlab="Clay %")

plot(data=horizons(mlra18.upper10mineral), estimated_om ~ sand,
     main='Organic Matter % - upper 10cm - MLRA 18', ylab="Organic Matter %",
     xlab="Sand %")
