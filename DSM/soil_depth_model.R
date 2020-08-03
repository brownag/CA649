# soil depth model
library(aqp)
library(soilDB)
library(raster)
library(magrittr)

f.all <- fetchNASIS()

# remove theresas TUDs with bad coordinates, and anything missing coordinates
f <- f.all
bad.idx <- which(is.na(f.all$x_std) | f.all$x_std > -119)
if(length(bad.idx))
  f <- f[-bad.idx,]

if(length(slot(f, 'sp')) == 1) {
  coordinates(f) <- ~ x_std + y_std
  proj4string(f) <- "+proj=longlat +datum=WGS84"
}

# inspect set spatially
plot(slot(f, 'sp'))

# inspect observed distribution in bedrock depth
plot(density(f$bedrckdepth, na.rm=T))

# NA bedrock depth in 4 cases
f$bedrckdepth[!is.finite(f$bedrckdepth)] <- NA
sum(is.na(f$bedrckdepth))

# basic model - predict bedrock depth using field-measured slope gradient
#   hypothesis is soils tend to be shallow to moderately deep
#    complex, steep, colluvially-active slopes slopes 
#    and material-accumulating areas have deeper soils
fm0 <- bedrckdepth ~ slope_field
m0 <- lm(data=site(f), formula = fm0)
summary(m0) # slope is a significant predictor at 0.05 alpha level

plot(data=site(f), fm0, ylim=c(0,200))
abline(m0)
conf <- predict(m0, newdata=data.frame(slope_field=0:90), interval = 'confidence')
lines(0:90, conf[,2], lty=2, col="RED")
lines(0:90, conf[,3], lty=2, col="RED")
pred <- predict(m0, newdata=data.frame(slope_field=0:90), interval = 'prediction')
lines(0:90, pred[,2], lty=2, col="BLUE")
lines(0:90, pred[,3], lty=2, col="blue")

# average depth on 0% slope is 53cm - confidence interval: 45 to 63cm
coef(m0)
# LOW ELEVATION MODEL
#  for each 10% increase in slope, there is ~3.5cm more soil material, on average
#   the confidence interval is 1.5 to 7cm
# HIGH ELEVATION MODEL
#  for each 10% increas in slope, there is ~5.5cm more soil material on averaage
#   the confidence interval is 3.5 to 8cm
confint(m0)

# define a function to break up the data into "slope phases"
split.points <- function(p, attr, lo, hi) {
  if(length(lo) != length(hi)) {
    stop('low and hi threshold vectors must be of same length')
  }
  do.call('cbind', lapply(as.list(1:length(lo)), function(i) {
    dplyr::between(p[[attr]], lo[i], hi[i])
  }))
}

# try dividing dataset int 0-15 and 15-90
slope.phases <- split.points(f, 'slope_field', c(0,15), c(15,90))

# same model, only just run on 0-15% phase
f0to15 <- site(f)[which(slope.phases[,1]),]
m0a <- lm(data=f0to15, fm0)
summary(m0a) 
# low elevation: slope is not a good predictor of soil depth for 0-15% slopes
# high elevation: same -- possibly deeper soils on lowest slopes as well

plot(data=f0to15, fm0, ylim=c(0,200))
abline(m0a)
conf <- predict(m0a, newdata=data.frame(slope_field=0:90), interval = 'confidence')
lines(0:90, conf[,2], lty=2, col="RED")
lines(0:90, conf[,3], lty=2, col="RED")
pred <- predict(m0a, newdata=data.frame(slope_field=0:90), interval = 'prediction')
lines(0:90, pred[,2], lty=2, col="BLUE")
lines(0:90, pred[,3], lty=2, col="blue")

# average depth on 0% slope is 59cm - confidence interval: 45 to 63cm
coef(m0a)
# LOW ELEVATION MODEL
#  for each 10% increase in slope, there is  more soil material, on average
#   the confidence interval is  to 
# HIGH ELEVATION MODEL 0 to 15
#  for each 10% increase in slope, there is ~2.5cm less soil material on averaage
#   the confidence interval is -15cm to +10cm (contains zero -- no relationship)
confint(m0a)

# same model, only just run on 15-90% phase
f15to90 <- site(f)[which(slope.phases[,2]),]
m0b <- lm(data=f15to90, fm0)
summary(m0b) # slope is a marginally significant predictor at 0.1 alpha level

plot(data=f15to90, fm0, ylim=c(0,200))
abline(m0b)

conf <- predict(m0b, newdata=data.frame(slope_field=0:90), interval = 'confidence')
lines(0:90, conf[,2], lty=2, col="RED")
lines(0:90, conf[,3], lty=2, col="RED")

pred <- predict(m0b, newdata=data.frame(slope_field=0:90), interval = 'prediction')
lines(0:90, pred[,2], lty=2, col="BLUE")
lines(0:90, pred[,3], lty=2, col="blue")

library(raster)
r <- stack(list.files(path = "Geodata/derived/", pattern = "tif", full.names = T))
sp.xtract <- spTransform(slot(f, 'sp'), CRS(proj4string(r)))
site(f) <- data.frame(peiid=profile_id(f), extract(r, sp.xtract))

site(f) <- profileApply(f, function(p) {
  
  # calculate depth to bedrock if present
  d <- estimateSoilDepth(p, p = "R|Cr|Cd|qm",
                         no.contact.depth = 0, 
                         no.contact.assigned = NA)
  
  # calculate maximum bottom depth described
  m <- max(p$hzdepb, na.rm=TRUE)
  
  # construct single row data frame result
  data.frame(peiid = profile_id(p), 
             bedrock_described = !is.na(d <= m), 
             soil_depth = min(d, m, na.rm=T))
  
}, frameify=TRUE)

plot(filter(f, !bedrock_described), label="pedon_id")

col.names <- c('peiid','x','y','soil_depth', 'depth.class', names(r))
predictors <- site(f)[complete.cases(site(f)[,col.names]), col.names]

m1 <- step(lm(soil_depth ~ . - peiid - aspect, data=predictors))
summary(m1)

plot(predict(m1, predictors), 
     predictors$soil_depth, 
     xlim = c(0, 200), ylim = c(0, 200))
abline(0,1)
abline(lm(predict(m1, predictors) ~ predictors$soil_depth))

survival.data <- site(f)[,c("soil_depth", 
                            "bedrock_described",
                            names(r))]

predictors$is_shallow <- FALSE
predictors$is_shallow[predictors$depth.class == "shallow"] <- TRUE
predictors$is_shallow[predictors$depth.class == "very.shallow"] <- TRUE
predictors$is_shallow <- factor(predictors$is_shallow)

rf0 <- randomForest(data = predictors, 
                    classwt = c(100000,1),
                    is_shallow ~ . - peiid - soil_depth - depth.class - x - y)

folds <- CAST::CreateSpacetimeFolds(predictors, spacevar = "peiid")

plot(as(extent(r), 'SpatialPolygons'))
i <- 1
lapply(folds$index, function(ff) {
  pedons <- aqp::filter(f, peiid %in% predictors$peiid[ff])
  points(spTransform(slot(pedons, 'sp'), CRS(proj4string(r))), pch=i)  
  i <<- i + 1
})

load("depth_ffs.Rda")
res <- CAST::ffs(predictors[,2:length(names(predictors))], response = predictors$soil_depth)
out <- res$finalModel

rtiny <- crop(r, extent(r) / 100)

foo <- raster::predict(rtiny, out)#, filename="ffs_prediction_depth_cm.tif")
plot(foo)

plot(predict(rf0)[as.numeric(names(m0$fitted.values))] ~  predict(m0), xlim = c(0, 200), ylim = c(0, 200))
m1 <- lm(predict(rf0)[as.numeric(names(m0$fitted.values))] ~ predict(m0))
abline(m1)
abline(0,1)

summary(m1)

plot(density(survival.data$soil_depth),
     col="blue", main="Observed versus Predicted (ordinary RF)", ylim=c(0,0.05))
lines(density(predict(rf0)))

abline(v=median(predict(rf0)), lwd=2)
abline(v=median(survival.data$soil_depth), col="blue", lwd=2)

abline(v=quantile(predict(rf0), probs=c(0.1, 0.9)), lty=2)
abline(v=quantile(survival.data$soil_depth, probs=c(0.1, 0.9)), lty=2, col="blue")

# lines(density(rnorm(1e6, mean(df$observed), sd(df$observed))))
# abline(v=quantile(df$observed), lty=2, col='blue')
# abline(v=quantile(pred[,1]), lty=2)
legend("topright", 
       c("Observed", "Predicted", 
         "Observed Median", "Predicted Median", 
         "Observed Q10/Q90", "Predicted Q10/Q90"), 
       col=c("BLUE","BLACK","BLUE","BLACK", "BLUE","BLACK"), 
       lty=c(1,1,1,1,2,2), lwd=c(1,1,2,2,1,1))

