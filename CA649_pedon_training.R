# CA649_pedon_training.R
library(aqp)
library(soilDB)
library(rgdal)

# get observations from CA630 + CA649
pedons <- fetchNASIS(SS = FALSE)
bad.idx <- which(is.na(pedons$x_std) | is.na(pedons$y_std))
if(length(bad.idx)) {
  pedons <- pedons[-bad.idx,]  
}

coordinates(pedons) <- ~ x_std + y_std
proj4string(pedons) <- '+proj=longlat +datum=WGS84'

ca630.legend <- soilDB::get_mapunit_from_SDA(WHERE = "areasymbol = 'CA630'")
ca630.legend.mvo <- ca630.legend[grepl(x = ca630.legend$muname, pattern = "^Loafercreek|^Bonanza|^Gopheridge|^Jasperpeak|^Gardellones|^Aquic"),]
ca630.legend.mvo <- ca630.legend.mvo[!grepl(x = ca630.legend.mvo$muname, pattern = "Urban"),]

chunk_SDA_spatial <- function(mukey.list, nchunk = 10) {
  mukey.chunk <- (1:length(mukey.list) %% nchunk) + 1
  s <- NULL
  
  for(i in 1:max(mukey.chunk)) {
    idx <- which(mukey.chunk == i)
    
    q <- paste0("SELECT G.MupolygonWktWgs84 as geom, mapunit.mukey FROM mapunit 
                CROSS APPLY SDA_Get_MupolygonWktWgs84_from_Mukey(mapunit.mukey) as G 
                WHERE mukey IN ", format_SQL_in_statement(mukey.list[idx]))
    
    #NB: FedData has a different (much simpler, but not equivalent) definition of SDA_query
    #    it also uses the post/rest interface
    sp.res.sub <- suppressMessages(soilDB::SDA_query(q))
    s.sub <- soilDB::processSDA_WKT(sp.res.sub)
    
    if(is.null(s)) {
      s <- s.sub
    } else {
      s <- rbind(s, s.sub)
    }
  }
  return(s)
}
ca630.spatial <- chunk_SDA_spatial(mukey.list = unique(ca630.legend.mvo$mukey))
ca630.spatial <- spTransform(ca630.spatial, CRSobj = CRS(proj4string(pedons)))
ca630.spatial <- merge(ca630.spatial, ca630.legend, all.x=TRUE, sort=F)
pedons$musym <- over(as(pedons, 'SpatialPointsDataFrame'), ca630.spatial)$musym
pedons$taxonname <- toupper(pedons$taxonname)
pedons <- pedons[which(!is.na(pedons$musym) | grepl(pedons$pedon_id, pattern = "CA043")),]

## inspect the data
groupedProfilePlot(pedons, groups='taxonname', color='total_frags_pct')
plot(density(pedons$total_frags_pct, na.rm=T))

###### ----------------------------------------
######  CREATE SITE-LEVEL GROUPING VARIABLE(S)
###### ----------------------------------------
pscs <- profileApply(pedons, estimatePSCS, simplify = FALSE)

df.pscs <- as.data.frame(do.call('rbind', pscs))
names(df.pscs) <- c("pscs_t","pscs_b")

depth.weighted.average <- function(spc, tdepth, bdepth, attr, ...) {
  #expand `attr` in formula
  custom.formula <- formula(paste0(tdepth,":",bdepth," ~ ", 
                                   paste0(attr, collapse=" + ")))
  # calculate a depth-weighted average using aqp::slice()
  return(mean(slice(spc, custom.formula, just.the.data=TRUE)[[attr]], na.rm = TRUE))
}

slot(pedons, 'site') <- cbind(site(pedons), df.pscs)

site(pedons)$pscs_clay <- profileApply(pedons, function(p) {
  return(depth.weighted.average(p, p$pscs_t, p$pscs_b, 'clay'))
})

site(pedons)$pscs_frags <- profileApply(pedons, function(p) {
  return(depth.weighted.average(p, p$pscs_t, p$pscs_b, 'total_frags_pct'))
})

site(pedons) <- getSoilDepthClass(pedons)

# create binary outcomes (not using these to _train_ the models at this point)
pedons$is_fine <- factor(pedons$pscs_clay >= 35)
pedons$is_skeletal <- factor(pedons$pscs_frags >= 35)
pedons$is_shallow <- factor(pedons$depth <= 50)

table(pedons$is_fine)
table(pedons$is_skeletal)
table(pedons$is_shallow)

plot(as(pedons,'SpatialPointsDataFrame'), col=viridis::viridis(2)[match(site(pedons)$is_shallow, c(T,F))])
