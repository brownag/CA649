# mvo model extended to CA649

# train random forest models to predict soil depth class, PSCS clay and PSCS fragments on CA630 extent
# these predictions (.tif) can be used with CA649_series_model.R to create maps of soil components

library(rgdal)
library(raster)
library(clhs)
library(rasterVis)

source('CA649_pedon_training.R')

# load covariates -- 10m CA630 data for training; 10m CA649 data for prediction
ca630_r <- stack('D:/stats2/dem_derivatives_mask.tif')
r <- stack('E:/CA649/ca649_terrain_derivatives_10m.tif')

covariate.names <- c("slope","aspect","cross_curv","long_curv","conv_idx",
                     "closed_depression","flow_acc","TWI","LS","channel_base",
                     "vdist_channel","valley_depth","relative_slope_pos",
                     "TRI","MRVBF","MRRTF","VRM")

### load training data
ca630_covariates <- as.data.frame(extract(ca630_r, slot(pedons, 'sp')))
colnames(ca630_covariates) <- covariate.names
names(ca630_r) <- covariate.names

# do conversions from degrees to percent slope
ca630_r$slope <- tan(ca630_r$slope) * 100
ca630_covariates$slope <- tan(ca630_covariates$slope) * 100

# add extracted data to SPC
site(pedons) <- cbind(site(pedons), ca630_covariates)

covariate.stack <- c('depth', 'pscs_clay', 'pscs_frags', colnames(ca630_covariates))
#pedons <- pedons[complete.cases(site(pedons)[,c('depth', 'pscs_clay', 'pscs_frags')]), ]

#pedons.clhs <- clhs(site(pedons)[,c('depth', 'pscs_clay', 'pscs_frags')], 
#                    iter=100, size=50, simple=TRUE)
pedons.clhs <- 1:length(pedons)
  
par(mfrow=c(1,1))
plot(pedons@sp[pedons.clhs])
points(pedons@sp)

plot(density(pedons$depth))
lines(density(pedons[pedons.clhs]$depth), col="BLUE")

plot(density(pedons$pscs_clay, na.rm=T))
lines(density(pedons[pedons.clhs]$pscs_clay, na.rm=T), col="BLUE")

plot(density(pedons$pscs_frags, na.rm=T))
lines(density(pedons[pedons.clhs]$pscs_frags, na.rm=T), col="BLUE")

# inspect extracted values
#  compare field measured versus raster extracted slope grade (%)
plot(pedons[pedons.clhs]$slope_field ~ pedons$slope[pedons.clhs], 
      xlim=c(0,90), ylim=c(0,90),
      xlab="Raster-extracted Slope, %", ylab="Field-measured Slope, %")
abline(0,1)
stack.resamp <- as.data.frame(sampleRegular(ca630_r, 10000))
clhs.samples.idx <- clhs(stack.resamp, iter=10, size=length(pedons))
clhs.samples <- stack.resamp[clhs.samples.idx,]

# create training and test dataframes
attr <- c('peiid', 'taxonname', 'depth', 
          'depth.class', 'pscs_clay', 'pscs_frags', 
          'is_shallow','is_fine', 'is_skeletal', 
          colnames(ca630_covariates))
train.idx <- sample(1:length(pedons[pedons.clhs]), 0.90*length(pedons[pedons.clhs]))
test.idx <- which(!(1:length(pedons[pedons.clhs]) %in% train.idx))
train <- site(pedons[pedons.clhs])[train.idx, attr]
test <- site(pedons[pedons.clhs])[test.idx, attr]
train.complete <- train[complete.cases(train),]

# fit a random forest model for each quantitative soil attribute
#   the models use different aspects of terrain shape and roughnes to predict:
#     1. soil depth
#     2. PSCS fragments
#     3. PSCS clay

# note that channel_base has been removed as it tends to be SSA-specific and produces odd results due
# to dependence on local topography and extent of input DEM (when terrain derivatives are made). 
# vdist_channel appears to have good explanatory power on its own and is less site-specific

depth_model <- randomForest::randomForest(depth ~ slope + conv_idx + closed_depression + 
                                            TWI + LS + vdist_channel + valley_depth + 
                                            relative_slope_pos + TRI + MRVBF + MRRTF + VRM, 
                                          data = train.complete, importance = TRUE, ntree=4000)

depth_class_model <- randomForest::randomForest(is_shallow ~ slope + conv_idx + closed_depression + 
                                            TWI + LS + vdist_channel + valley_depth + 
                                            relative_slope_pos + TRI + MRVBF + MRRTF + VRM, 
                                          data = train.complete, importance = TRUE, ntree=4000,
                                          classwt=c(0,2), strata=train.complete$is_shallow)

library(rpart.plot)
library(caret)
library(ggplot2)
library(AppliedPredictiveModeling)

V <- 10
T <- 4
Seed <- 1

TrControl <- trainControl(method = "repeatedcv",
                          number = V,
                          repeats = T)

frag_class_model <- randomForest::randomForest(is_skeletal ~ slope + conv_idx + closed_depression + 
                                           TWI + LS + vdist_channel + valley_depth + 
                                           relative_slope_pos + TRI + MRVBF + MRRTF + VRM, 
                                         data = train.complete, importance = TRUE, ntree=4000,
                                         classwt=c(0,2), strata=train.complete$is_skeletal)

clay_model <- randomForest::randomForest(pscs_clay ~ slope + conv_idx + closed_depression + 
                                           TWI + LS + vdist_channel + valley_depth + 
                                           relative_slope_pos + TRI + MRVBF + MRRTF + VRM, 
                                         data = train.complete, importance = TRUE, ntree=4000)

clay_class_model <- randomForest::randomForest(is_fine ~ slope + conv_idx + closed_depression + 
                                           TWI + LS + vdist_channel + valley_depth + 
                                           relative_slope_pos + TRI + MRVBF + MRRTF + VRM, 
                                         data = train.complete, importance = TRUE, ntree=4000,
                                         classwt=c(0.01,2), strata=train.complete$is_fine)

# here we save the randomForest model objects, the training data, and the clhs samples (for comparison)
save(file = 'CA630_rf_models.Rda', 
     list = c('depth_model','frag_model','clay_model', 
              'depth_class_model','frag_class_model','clay_class_model', 
              'train','test','train.complete', 
              'clhs.samples'))
# in the future, when trying to make predictions from 

load(file = 'CA630_rf_models.Rda')

# ###### --------------------------------------------------
# ######  VARIABLE IMPORTANCE
# ###### --------------------------------------------------
randomForest::varImpPlot(depth_model)
randomForest::varImpPlot(depth_class_model)

randomForest::varImpPlot(frag_model)
randomForest::varImpPlot(frag_class_model)

randomForest::varImpPlot(clay_model)
randomForest::varImpPlot(clay_class_model)

# ###### --------------------------------------------------
# ######  MAKING PREDICTIONS USING TEST DATA (and compare to training data)
# ###### --------------------------------------------------

#par(mfrow=c(1,3))
# depth class model
depth_train <- predict(depth_model, newdata = train)
depth_test <- predict(depth_model, newdata = test)
plot(train$depth, depth_train, 
     xlim=c(0,200),ylim=c(0,200),
     xlab="Field Depth, cm", ylab="Predicted Depth, cm")
points(test$depth, depth_test, pch=3)
abline(0,1)
legend(0, 200, legend = c("Training points", "Test points"), pch=c(1,3))

# pscs rock fragment model
frag_train <- predict(frag_model, newdata = train)
frag_test <- predict(frag_model, newdata = test)
plot(train$pscs_frags, frag_train, 
     xlim=c(0,80),ylim=c(0,80),
     xlab="Field PSCS Fragments, %", ylab="Predicted PSCS Fragments, %")
points(test$pscs_frags, frag_test, pch=3)
abline(0,1)
legend(0, 80, legend = c("Training points", "Test points"), pch=c(1,3))

# pscs clay model
clay_train <- predict(clay_model, newdata = train)
clay_test <- predict(clay_model, newdata = test)
plot(train$pscs_clay, clay_train, 
     xlim=c(18,60),ylim=c(18,60),
     xlab="Field PSCS Clay, %", ylab="Predicted PSCS Clay, %")
points(test$pscs_clay, clay_test, pch=3)
abline(0,1)
legend(18, 60, legend = c("Training points", "Test points"), pch=c(1,3))

names(r) <- covariate.names

### look at classification accuracy
depth_class_model
frag_class_model
clay_class_model

# convert to %
r$slope <- tan(r$slope) * 100

depth_pred <- predict(r, depth_model, progress="text")
p <- levelplot(depth_pred, main="Predicted Depth to Contact (cm)")
plot(p)
writeRaster(depth_pred,'CA649_predicted_depth_cm.tif')

frag_pred <- predict(r, frag_model, progress="text")
p <- levelplot(frag_pred, main="Predicted PSCS Fragments (%)")
plot(p)
writeRaster(frag_pred,'CA649_predicted_frag_pct.tif')

clay_pred <- predict(r, clay_model, progress="text")
p <- levelplot(clay_pred, main="Predicted PSCS Clay (%)")
plot(p)
writeRaster(clay_pred,'CA649_predicted_clay_pct.tif')
