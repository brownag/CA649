# CA649 series model

# this uses rasters of predicted soil depth, PSCS fragments and PSCS clay to assign series name or concept

library(raster)

depth_pred <- raster('CA649_predicted_depth_cm.tif')
frag_pred<- raster('CA649_predicted_frag_pct.tif')
clay_pred <- raster('CA649_predicted_clay_pct.tif')

pred <- stack(depth_pred,frag_pred,clay_pred)
names(pred) <- c("depth","pscs_frag","pscs_clay")

series.names <- c("Bonanza","Jasperpeak","Clayey, Shallow","CL-SK, Shallow",
                  "Loafercreek","Gopheridge","Argonaut", "Cl-SK, ModDeep",
                  "Motherlode","Gardellones","Fine, Deep", "Cl-Sk, Deep")
series.id <- 1:length(series.names)
names(series.id) <- series.names

series_fun <- function(depth,pscs_frag,pscs_clay) {
  shallow.idx <- depth <= 50
  moddeep.idx <- depth > 50 & depth <= 100
  deep.idx <- depth >= 100
  skeletal.idx <- !(pscs_frag < 35)
  fine.idx <- !(pscs_clay < 35)
  nonfine <- !fine.idx
  buf <- rep(NA, length(depth))
  buf[shallow.idx & !skeletal.idx & !fine.idx] <- 1
  buf[shallow.idx & skeletal.idx & !fine.idx] <- 2
  buf[shallow.idx & !skeletal.idx & fine.idx] <- 3
  buf[shallow.idx & skeletal.idx & fine.idx] <- 4
  buf[moddeep.idx & !skeletal.idx & !fine.idx] <- 5
  buf[moddeep.idx & skeletal.idx & !fine.idx] <- 6
  buf[moddeep.idx & !skeletal.idx & fine.idx] <- 7
  buf[moddeep.idx & skeletal.idx & fine.idx] <- 8
  buf[deep.idx & !skeletal.idx & !fine.idx] <- 9
  buf[deep.idx & skeletal.idx & !fine.idx] <- 10  
  buf[deep.idx & !skeletal.idx & fine.idx] <- 11
  buf[deep.idx & skeletal.idx & fine.idx] <- 12
  return(buf)
}

foo <- overlay(pred, unstack=TRUE, fun = series_fun, 
               filename = "CA649_series_prediction.tif",
               progress = "text")
