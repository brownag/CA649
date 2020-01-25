#CA649 ecosite model
library(aqp)
library(soilDB)
library(raster)
library(sp)

source("config.R")

ca630.legend <- soilDB::get_mapunit_from_SDA(WHERE = "areasymbol = 'CA630'")
ca649.legend <- soilDB::get_mapunit_from_SDA(WHERE = "areasymbol = 'CA649'")

unique.mukey <- unique(c(ca630.legend$mukey, ca649.legend$mukey))

all.comp <- soilDB::get_component_from_SDA(WHERE = paste0('mukey IN ',
                                              format_SQL_in_statement(unique.mukey)),
                                           duplicates = TRUE)

# update ecosite assignments from daves soil sort
mlra18.mukey.to.pESD <- read.csv("Ecosites/MLRA18_coiid_to_pESD.csv")

all.comp <- merge(all.comp, mlra18.mukey.to.pESD[,-1], by="mukey", all.x=T, sort=F)

train.map <- do.call('rbind',
               lapply(split(all.comp, f=all.comp$mukey),
 function(dmu) {
  with(dmu, {
    eco.comp <- aggregate(comppct_r, by=list(ecoclassidPES), FUN=sum)
    res <- eco.comp[order(eco.comp[,2], decreasing = TRUE),]
    res[1,]
  })
 }))
train.map$mukey <- rownames(train.map)
rownames(train.map) <- NULL
colnames(train.map) <- c("ecoclassid","comppct_total","mukey")

# get 18XI ecosite dominant mapunits--excluding 999
res.18xi <- merge(train.map, all.comp[,c("mukey","nationalmusym")], all.x=TRUE, sort=FALSE)
res.18xi <- res.18xi[!grepl(res.18xi$ecoclassid, pattern = "999"),]

training.mus <- unique(res.18xi[!is.na(res.18xi$ecoclassid),])

# 
# mu <- SDA_spatialChunk(training.mukey)
# mu.set <- unique(training.mukey)
# slot(mu, 'data') <- merge(slot(mu, 'data'), training.mus, sort=F, all.x=T)
# col.lut <- viridis::viridis(length(unique(mu$ecoclassid)))
# names(col.lut) <- unique(mu$ecoclassid)
# newnames <- names(mu)
# newnames[which(newnames=='gid')] <- "pID"
# names(mu) <- newnames
# # visualize extent of different dominant ecosites (>=65% aggregate)
# plot(mu, col=col.lut[match(s$ecoclassid, names(col.lut))], border=NA)
# 
# # convert to projected 
# mu <- sp::spTransform(mu, sp::CRS("+proj=aea +lat_1=30 +lat_2=50 +lat_0=40 +lon_0=-125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
# 
# # write shapefile
# mu$comppct_total <- as.numeric(mu$comppct_total)
# rgdal::writeOGR(mu, dsn='.', 'dominant_ecosite_XI', 'ESRI Shapefile')


mu <- rgdal::readOGR(dsn='.', 'dominant_ecosite_XI')
mu <- sp::spTransform(mu, sp::CRS("+proj=aea +lat_1=30 +lat_2=50 +lat_0=40 +lon_0=-125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

names(mu) <- c('mukey','gid','nmusym','ecosite','comppct')
mu.set <- unique(mu$nmusym)
mu$pID <- mu$gid

# print("Sampling raster stack by mapunit...")
# .timer.start <- Sys.time()
# sampling.res <- suppressWarnings(sharpshootR::sampleRasterStackByMU(mu, mu.set, 'nmusym', raster.list, pts.per.acre = pts.per.acre, estimateEffectiveSampleSize = FALSE))
# .timer.stop <- Sys.time()
# .sampling.time <- format(difftime(.timer.stop, .timer.start, units='mins'), digits=2)
# print(paste0("Completed sampling of raster stack for poly symbols (",paste(mu.set,collapse=","),") in ",.sampling.time,"."))
# 
# sampling.res$raster.samples$.id <- factor(sampling.res$raster.samples$.id, levels=mu.set)
# print("DONE!")
# 
# 
# foo <- sampling.res$raster.samples
# res3 <- stats::reshape(foo, timevar=c("variable","variable.type"),
#                v.names = c("value"),
#                idvar = c(".id","pID","sid"),
#                direction="wide")
# names(res) <- c(".id", "variable.type", "pID", "sid", "maatC", 
#                 "mapmm", "effpmm", 
#                 "ffd", "gdd", 
#                 "elev", "slope", "abr", 
#                 "twi", "geomorphon", "curvcl", 
#                 "nlcd", "nlcd_crop", 
#                 "aspect")
# save(res, file = "ecosite_samples.Rda")

load(file = "ecosite_samples.Rda")

# cook down samples with cLHS
clhs.sample.idx <- clhs::clhs(res[,!names(res) %in% c(".id","variable.type","pID","sid")], size=round(0.01*nrow(res)), iter=1000)

clhs.samples <- res[clhs.sample.idx, ]

nmusym.to.eco.lut <- as.character(training.mus$ecoclassid)
names(nmusym.to.eco.lut) <- as.character(training.mus$nationalmusym)

clhs.samples$ecoclassid <- nmusym.to.eco.lut[as.character(clhs.samples$.id)]
training.data <- clhs.samples
training.data <- training.data[training.data$ecoclassid %in% c("F018XC201CA", "F018XI200CA", "F018XI201CA", "F018XI202CA", "R018XC105CA", "R018XI106CA",  "R018XI163CA"),]
training.data$ecoclassid <- factor(as.character(training.data$ecoclassid))
training.data$mlraid <- factor(substr(training.data$ecoclassid, 2, 4), 
                                 levels=c("017","018","22A"))
training.data$lruid <- factor(substr(training.data$ecoclassid, 5, 6), 
                               levels=c("XI","XC"))
training.data$ecogrpid <- factor(substr(training.data$ecoclassid, 1, 1), 
                                 levels=c("R","F"))
training.data$ecoclassid <- factor(training.data$ecoclassid)

save(training.data, file = "ecosite_training_data.Rda")
load(file = "ecosite_training_data.Rda")

training.data <- training.data[complete.cases(training.data),]

# just use continuous predictors for rf1
rf1 <- randomForest::randomForest(ecoclassid ~ maatC+mapmm+effpmm+abr+slope, data=training.data)

save(rf1, file = "ecoclass_continuous_rf.Rda")
load(file = "ecoclass_continuous_rf.Rda")

pred.pt <- predict(rf1, newdata = training.data)
head(pred.pt)

# now build rasterstack for spatial prediction
alignRasters <- function(rasters, categorical, project.extent) {
  target <- rasters[[1]]
  target <- projectRaster(from = target, crs=proj4string(project.extent))
  target <- crop(target, extent(project.extent))
  message("converted template raster to project extent CRS and cropped")
  res <- list(target)
  for(i in 2:length(rasters)) {
    tmp <- crop(rasters[[i]], spTransform(project.extent, CRS(proj4string(rasters[[i]]))))
    if(!categorical[i])
      tmp <- raster::projectRaster(from = rasters[[i]], to = target, method = "bilinear")
    else
      tmp <- raster::projectRaster(from = rasters[[i]], to = target, method = "ngb")
    tmp <- crop(tmp, extent(target))
    res[[i]] <- tmp
    message(paste0("Aligned ", names(rasters[[i]]), " with ", names(target),"."))
  }
  return(res)
}
ca649_b <- rgdal::readOGR("E:/CA649/Geodata/source/CA649_b.shp")
input.ras <- lapply(as.list(unlist(raster.list)), raster)
predictor.stack <- alignRasters(input.ras, 
                                c(F,F,F,F,F,F,F,F,F,T,T,T,T,F), 
                                project.extent = ca649_b)
predictor.stack <- stack(predictor.stack)
names(predictor.stack) <- c("maatC", 
                "mapmm", "effpmm",
                "ffd", "gdd",
                "elev", "slope", "abr",
                "twi", "geomorphon", "curvcl",
                "nlcd", "nlcd_crop",
                "aspect")
writeRaster(predictor.stack, filename = "ca649_predictor_stack_800m.tif")

prediction1 <- raster::predict(model=rf1, object=predictor.stack, 
                               filename="prediction1.tif", 
                               overwrite=TRUE, progress="text")
plot(prediction1)
plot(ca649_b, add=T)

lut <- data.frame(id=1:length(rf1$classes), ecositeclass=(rf1$classes))
lut
