# predict mu suites
library(rgdal)
library(raster)
library(fasterize)

training.shp <- readOGR("../CA649/Premap/ecoboundaries_v1.shp")
predictor.stack <- stack('predictor_stack_800m.tif')
names(predictor.stack) <- c("maatC", "mapmm", "effpmm",
                            "ffd", "gdd",
                            "elev", "slope", "abr",
                            "twi", "geomorphon", "curvcl",
                            "nlcd", "nlcd_crop")

training.shp <- spTransform(training.shp, CRS(proj4string(predictor.stack)))
training.mask <- fasterize::fasterize(sf::st_as_sf(training.shp), predictor.stack[[1]])
predictor.stack <- mask(predictor.stack, training.mask)

plot(predictor.stack)

all.samples <- sampleRegular(predictor.stack, size = 100000, sp=TRUE)
all.samples <- all.samples[complete.cases(slot(all.samples,'data')),]


clhs.idx <- clhs::clhs(all.samples[,1:9], size=round(0.1*nrow(all.samples)))
clhs.idx$class <- over(clhs.idx, training.shp)$Name
clhs.idx$class <- factor(clhs.idx$class)#, levels=c("LowElevFoothills","SouthernFoothills",
                                       #                  "IntermediateFoothills","HighHillsRidges",

complete.idx <- which(complete.cases(slot(clhs.idx, 'data')))
clhs.idx <- clhs.idx[complete.idx,]
train.idx <- sample(1:nrow(clhs.idx), round(0.9*nrow(clhs.idx)))
train.samples <- clhs.idx[train.idx,]
test.samples <- clhs.idx[-train.idx,]

rf <- randomForest::randomForest(class ~ slope + abr + twi + effpmm + ffd, data=slot(train.samples, 'data'))

prediction2 <- raster::predict(model=rf, object=stack(predictor.stack),
                               filename="prediction_musuite.tif", overwrite=T)
spplot(prediction2)
rf$classes
