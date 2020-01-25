### MU GIS Summary Report
### 2018-05-29
### D.E. Beaudette and J. Wood
###
### configuration file, edit as needed
###

#########################
### Raster Data Sources #
#########################

# data sources can be "commented-out" using the "#" character
# be sure that there is no trailing "," after the last item in each list
# raster summaries are displayed in the same order in which they are listed below

raster.list <- list(
  continuous=list(
    `Mean Annual Air Temperature (degrees C)`='L:/NRCS/MLRAShared/Geodata/project_data/MUSum_PRISM/final_MAAT_800m.tif',
    `Mean Annual Precipitation (mm)`='L:/NRCS/MLRAShared/Geodata/project_data/MUSum_PRISM/final_MAP_mm_800m.tif',
    #`Moisture Deficit Anomaly (1981-2015)`='snf_mdeficit_anomaly_1981-2015.tif',
    `Effective Precipitation (mm)`='L:/NRCS/MLRAShared/Geodata/project_data/MUSum_PRISM/effective_precipitation_800m.tif',
    `Frost-Free Days`='L:/NRCS/MLRAShared/Geodata/project_data/MUSum_PRISM/ffd_50_pct_800m.tif',
    `Growing Degree Days (degrees C)`='L:/NRCS/MLRAShared/Geodata/project_data/MUSum_PRISM/gdd_mean_800m.tif',
    `Elevation (m)`='L:/NRCS/MLRAShared/Geodata/project_data/MUSum_30m_SSR2/DEM_30m_SSR2.tif',
    `Slope Gradient (%)`='L:/NRCS/MLRAShared/Geodata/project_data/MUSum_30m_SSR2/Slope_30m_SSR2.tif',
    #`NDVI (July 2018)`='LC08_L1TP_042034_20180705_20180717_01_T1_sr_ndvi_utm.tif',
    `Annual Beam Radiance (MJ/sq.m)`='L:/NRCS/MLRAShared/Geodata/project_data/ssro2_ann_beam_rad_int.tif',
    #`Compound Topographic Index`='L:/NRCS/MLRAShared/Geodata/project_data/ssro2_tci_int.tif',
    `SAGA TWI`='L:/NRCS/MLRAShared/Geodata/project_data/ssro2_saga_twi_int.tif'
    #`NLCD Impervious Surface (%))`='L:/NRCS/MLRAShared/Geodata/project_data/nlcd_impervious_2011_cropped.tif'
  ),
  categorical=list(
    `Geomorphon Landforms`='L:/NRCS/MLRAShared/Geodata/project_data/MUSum_Geomorphon/forms30_region2.tif',
    `Curvature Classes`='L:/NRCS/MLRAShared/Geodata/project_data/MUSum_Curvature/curvature_classes_30_class_region2.tif',
    `NLCD (2011)`='L:/NRCS/MLRAShared/Geodata/project_data/nlcd_2011_cropped.tif',
    `NASS Cropland Data Layer (2017)`='L:/NRCS/MLRAShared/Geodata/project_data/NASS/Ver2017_30m_cdls_clip.img'
  ),
  circular=list(
    `Slope Aspect (degrees)`='L:/NRCS/MLRAShared/Geodata/project_data/MUSum_30m_SSR2/Aspect_30m_SSR2.tif'
  )
)

pts.per.acre <- 2

# the most important quantiles (percentiles / 100) are: 0.1, 0.5 (median), and 0.9
# optionally reduce the number of quantiles for narrower tables
p.quantiles <- c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)


scaleDensityCurves <- TRUE


#####################################################################################
### output file names (OPTIONAL; uncomment to override defaults)                    #
### default will include a file-specific prefix and full list of MUSYMs summarized  #
### do not include .shp extension for shapefiles; automatically added by writeOGR() #
#####################################################################################
# shp.unsampled.fname <- 'un-sampled-polygons'# shapefile containing any unsampled polygons (usually too small or odd shape)
# shp.stats.fname <- 'polygons-with-stats' # shapefile containing median values / most likely classes by delineation
# shp.qc.fname <- 'poly-qc' # shapefile containing "proportion of samples outside 5-95% quantile range" by delineation

# csv.stats.fname <- 'poly-stats.csv' # comma-separated value file containing median values / most likely classes by delineation
# csv.qc.fname <- 'poly-qc.csv' # comma-separated value file containing "proportion of samples outside 5-95% quantile range" by delineation

########################################################
### Add estimate of confidence to box and whisker plots ###
########################################################

#enabling this feature will double the run time 
#enabling this feature will add "notches" to box and whisker plots
# that are close approximations to a confidence interval around the median
# adjusted for spatial autocorrelation
correct.sample.size <- FALSE


###########################################
### save samples after report has run ? ###
###########################################

# used for tinkering with a report .Rmd and debugging
# this will save samples to a file and subsequent report runs will use the saved samples
# not recommended for routine operation
cache.samples <- FALSE

