# dem-to-block_diagram.R
# 03/11/2019
# last revised: 04/05/2019
# @author: andrew brown; 
#          based on demo by dylan beaudette

library(rayshader)
library(rgl)
library(raster)
library(rgdal)
library(FedData)
library(viridis)
library(imager)
library(sf)
library(fasterize)

######
### SETUP
######

## 1. read shapefile for overlay (must cover full extent of elevation .TIF)
#       for example, ssurgo data symbolized on musym
#thematic_shp <- readOGR(dsn='L:/NRCS/MLRAShared/CA630/Archived_OFFICIAL_DB/final/FG_CA630_GIS_2018.gdb',
#                        layer="ca630_a")
# alternately, read in thematic raster directly
theme <- fasterize::raster('E:/CA649/LaPalomaHRC_cLHS_membership.tif')

## 2. thematic attribute - the column name in shapefile attribute table 
mu.col <- "MUSYM"

## 3. digital elevation model (TIFF, or other raster-compatible format) for a chunk of space
#     e.g. pan to desired area in ArcMap, and Data > Export Data > By Data Frame
elev <- fasterize::raster('S:/NRCS/Archive_Andrew_Brown/CA649/exemplars/LaPalomaHRC/LaPalomaHRC_dem.tif')
elev <- crop(elev,theme)
#elev <- raster('elev.tif')
theme <- resample(theme, elev)

## 4. OPTIONAL: resample raster input
##      if needed, you can resample to some other resolution/grid size
##      note that the DEM/hillshade and any derived overlays/shadows will be in this resolution

## copy elevation raster
# elev_template <- elev

## change raster resolution in template (leaving all else the same)
# res(elev_template) <- c(5, 5) #5 m by 5 m

## resample elevation raster to desired template 
# elev <- resample(elev, elev_template)

## save to file
# writeRaster(elev, filename='lidar_Tm_test_5m.tif')

######
### END SETUP
######

# convert elevation raster -> matrix
elmat <- matrix(extract(elev, extent(elev), buffer=1000), nrow=ncol(elev), ncol=nrow(elev))

# calculate (rectangular) boundary of DEM, use that to cut the overlay shapefile

# note: there is no specific reason your extent polygon has to be rectangular
#       but it is done here because rasters are commonly rectangular and we wan 
# extent.poly <- FedData::polygon_from_extent(elev)
# extent.poly <- spTransform(extent.poly, CRS(proj4string(thematic_shp)))
# thematic_shp <- crop(thematic_shp, y = extent.poly)
# 
# # assign numeric value that is 1:1 with mukey
# thematic_shp$munum <- match(thematic_shp[[mu.col]], unique(thematic_shp[[mu.col]]))

# produce raster
#thematic_shp <- spTransform(thematic_shp, CRS(proj4string(elev)))
#theme <- rasterize(x = thematic_shp, y = elev, 'munum')
#theme <- fasterize(sf = st_as_sf(thematic_shp), raster = elev, field = 'munum')

# inspect raster representation of theme musym
plot(theme)

# create RGB array from rasterized theme
tf <- "temp.png"
sgdf <- as(theme, 'SpatialGridDataFrame')
scl <- function(x) (x - min(na.omit(x))) / diff(range(na.omit(x))) 

# OPTIONAL: set SPECIFIC COLORS for certain theme levels
#n.grp <- length(unique(values(theme)))
#first.colors <- viridis(n.grp)
#first.colors[5] <- "#000000" # set 5th symbol to black
## ETC.
#alt.cols <- col2rgb(first.colors[values(theme)])
#cols <- alt.cols

# select some colors that ideally span the color ramps, get RGB
cols <- col2rgb(rev(viridis(256))[scl(values(theme)) * 255 + 1])
# comment this line out if you are setting colors manually with alt.colors

# get the RGB channels and put them in the SGDF
sgdf$red <- cols[1,]
sgdf$grn <- cols[2,]
sgdf$blu <- cols[3,]

# write to PNG then read it back with imager... to get the color array
# this is convenient but theoretically just rescaling the RGB channels in the SGDF
# would be sufficient (i.e. no file output required per se)
writeGDAL(sgdf[c("red", "grn", "blu")], tf, type="Byte", mvFlag=255, drivername="PNG")
load.array <- as.array(load.image(tf))

# take just the R, G and B -- and shift around the dimensions a bit for rayshader
my.array <- load.array[,,1,1:3]

# compute shadows
raymat <- ray_shade(elmat)
ambmat <- ambient_shade(elmat)

# with big files save the intermediates in case rgl crashes R or something
#save(raymat, ambmat, file = "intermediates.Rda")

# example: add overlay to static map
# elmat %>%
#   sphere_shade(texture = "desert") %>%
#   add_overlay(my.array, alphacolor=1) %>%
#   plot_map()

# INTERACTIVE 3D PLOT WITH RGL

# set perspective with right-mouse + drag
# zoom with mouse wheel
# rotate with left-mouse + drag

# interactive 3D plot via rgl
elmat %>%
  sphere_shade(texture = "bw") %>%
  add_overlay(my.array, alphacolor=1) %>%
  #add_water(detect_water(elmat), color="desert") %>%
  add_shadow(raymat, max_darken = 0.4) %>%
  add_shadow(ambmat, max_darken = 0.4) %>%
  plot_3d(elmat, zscale=4, fov=0, theta=30, water = 0,
          zoom=0.75, phi=45, windowsize = c(1000,800), lineantialias = TRUE)

# take a static picture of the rgl window
render_snapshot()

# important to clear the previous rgl window if any settings are adjusted
rgl::rgl.clear()
