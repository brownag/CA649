# build rasterstacks at various resolutions for MLRA-wide assessment

library(sp)
library(rgdal)
library(raster)

source('config.R')

# zoning of MLRA by ecosite groups/complexes; defines the target CRS
project.shp <- readOGR("../CA649/Premap/ecoboundaries_v1.shp")

# the target raster attributes (resolution) are defined by first element of raster.list

# # build rasterstack for spatial prediction
alignRasters <- function(rasters, categorical, project.extent) {
  
  # template raster is first one in list
  target <- rasters[[1]]
  
  # convert to CRS of project extent and crop
  target <- projectRaster(from = target, crs=proj4string(project.extent))
  target <- crop(target, extent(project.extent))
  message("converted template raster to project extent CRS and cropped")
  
  # loop through rasters 2:n and crop/reproject to project.extent
  res <- list(target)
  for(i in 2:length(rasters)) {
    tmp <- crop(rasters[[i]], spTransform(project.extent, CRS(proj4string(rasters[[i]]))))
    
    if(!categorical[i]) {
      # continuous rasters are interpolated by bilinear method
      tmp <- raster::projectRaster(from = rasters[[i]], to = target, method = "bilinear")
    } else {
      # categorical rasters are interpolated by nearest neighbor method
      tmp <- raster::projectRaster(from = rasters[[i]], to = target, method = "ngb")
    }
    
    # final crop of reprojected raster to match target raster
    tmp <- crop(tmp, extent(target))
    res[[i]] <- tmp
    message(paste0("Aligned ", names(rasters[[i]]), " with ", names(target),"."))
  }
  
  #return result (a "stackable" list of RasterLayers)
  return(res)
}

input.ras <- lapply(as.list(unlist(raster.list)), raster)
predictor.stack <- alignRasters(input.ras,
                                c(F,F,F,F,F,F,F,F,F,T,T,T,T,F),
                                project.extent = project.shp)
predictor.stack <- stack(predictor.stack)

writeRaster(predictor.stack, filename = "predictor_stack_800m.tif")
# names(predictor.stack) <- c("maatC",
#                             "mapmm", "effpmm",
#                             "ffd", "gdd",
#                             "elev", "slope", "abr",
#                             "twi", "geomorphon", "curvcl",
#                             "nlcd", "nlcd_crop")
